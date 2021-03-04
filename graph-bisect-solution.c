/** Print the cutsize and the imbalance of a graph partition
 * @param N the number of vertices in the graph
 * @param adjBeg beginning pointer of the adjacency list of each node
 * @param adj the list of adjacent vertex indices
 * @param part the array indicating the part index of each vertex (0 or 1)
 */
void printCutSizeImbal(
    int N,
    int *adjBeg,
    int *adj,
    int *part);

/* Return the futur label of the vertex
* @param vertex we compute its label
*/
int futurLabel(int vertex, int *adjBeg, int *adj, int *part){
  int Partedge[2] = {0,0}; // the two counters
  for (int k = adjBeg[vertex]; k < adjBeg[vertex+1]; k++){
    if (part[adj[k]]==0){
      Partedge[0]++;
    }else{
      Partedge[1]++;
    }
  } 
  if (Partedge[1-part[vertex]] > Partedge[part[vertex]] ){
    return 1 - part[vertex]; // strictly more in the other part means a change
  } else {
    return part[vertex]; // we keep the same label
  }
}

/*Return the rank of the processor in charge of the vertex */
int vertexToProc(int vertex, int N, int numProcs){
  return vertex * numProcs / N;
}

/* Return +1 for 1 and -1 for 0*/
int convert(int label){
  return 2*label -1;
}

/*Put the index of the beginnig of parts thanks to the length of them
* @param arrayProc a array of size numProcs
* @param arrayCumProc array of size numProcs+1 with the cumulative sum of arrayProc*/
void cumulativeSum(int *arrayProc, int *arrayCumProc, int numProcs){
  arrayCumProc[0]=0;
  for (int i = 1; i<numProcs+1; i++){
    arrayCumProc[i] = arrayCumProc[i-1] + arrayProc[i-1];
  }
}

/* Put in lenPart the number of element with label 0 or 1
* @param lenPart : array of size two lenPart[0 (resp 1)] = number of label 0 (resp 1)
*/
void repartionPart(int N, int *part, int *lenPart){
  lenPart[0]=0;
  lenPart[1]=0;
  for (int i = 0; i <N ; i++ ){
      if (part[i]==0){
        lenPart[0]++;
      } else if (part[i] == 1){
        lenPart[1]++;
      }
    }
}


int calculDeltaPart(int N,int* part){
  int delta = 0;
  for (int i = 0; i<N; i++){
    delta += (2*part[i]-1);
  }
  return delta;
}

/* Update recvCount thanks to the graph and an array vertexNeeded
* @param vertexNeeded : array of "bool" if the vertex is a neighboor of a vertex of the currentproc
*/
void detectNeighbour(int N, int *adjBeg, int *adj, int vertex, int *vertexNeeded, int *recvCount, int numProcs){
  int currentVertex;
  for (int k = adjBeg[vertex]; k < adjBeg[vertex+1]; k++){
    currentVertex = adj[k];
    if (vertexNeeded[currentVertex]==0){ // first need of this vertex
      vertexNeeded[currentVertex] = 1;
      recvCount[vertexToProc(currentVertex, N, numProcs)]++; //increase the number of vertices received from the proc in charge of this vertex
    }
  }

}

/* Compute recvIdx array for a given vertex
@param recvIdxBeg similar to adjBeg for recvIdx
*/

void assignNeighbour(int N,int *adjBeg,int *adj,int v,int* vertexNeeded, int *currentIdx, int *recvIdx, int *recvIdxBeg, int numProcs){
  int currentVertex;
  for (int k = adjBeg[v]; k<adjBeg[v+1]; k++){
    currentVertex = adj[k];
    if (vertexNeeded[currentVertex]==0){
      vertexNeeded[currentVertex] = 1;
      int procK = vertexToProc(currentVertex,N, numProcs);
      recvIdx[recvIdxBeg[procK]+ currentIdx[procK]] = currentVertex;
      currentIdx[procK] ++;
    }
  }
}

/** Computing the cutsize of a graph partition
 * @param N the number of vertices in the graph
 * @param adjBeg beginning pointer of the adjacency list of each node
 * @param adj the list of adjacent vertex indices
 * @param part the array indicating the part index of each vertex (0 or 1)
 * @param maxIters the number of iterations for which the algorithm must run
 * @param epsilon the load imbalance parameter >= 1.0
 * @param algoName the name of the algorithm to be executed
 */
void bisectGraph(
    int N,
    int *adjBeg,
    int *adj,
    int *part,
    int maxIters,
    double epsilon,
    char *algoName)
{
  int procRank, numProcs;
  printCutSizeImbal(N, adjBeg, adj, part);
  MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);


  if (strcmp(algoName, "bisect-seq") == 0) { // Sequential graph bisection
    if (procRank == 0) {
      if (numProcs != 1) {
        printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", numProcs);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    int otherPart;
    int lenPart[2]={0, 0}; // to compute the size of each label
    repartionPart(N, part,lenPart); // initialize lenPart
    float sizePartMax = (N/2)*epsilon;

    int *nextPart; 
    nextPart = (int *) malloc(N * sizeof(part[0]));

    for (int i = 1; i<= maxIters; i++){
      for (int j = 0; j< N; j++){
        otherPart = 1 - part[j];
        if (lenPart[otherPart] > sizePartMax){
          nextPart[j] = part[j];
        } else{
          nextPart[j] = futurLabel(j, adjBeg,adj, part );
          if (nextPart[j] != part[j]){
            lenPart[part[j]]--;
            lenPart[nextPart[j]]++;
          }
        }

      }
      memcpy(part, nextPart, sizeof(int)*N);

      printCutSizeImbal(N, adjBeg, adj, part);
    }
    free(nextPart);


  } else if (strcmp(algoName, "bisect-a2a") == 0) { // Graph bisection with all-to-all communication

      int signPart;
      int deltaTot = 0;
      int nextDelta = 0;
      int deltaLocal = 0;

      int lenVertex = N /numProcs;
      int offsetProc = lenVertex*procRank;
      int currentVertex;
      int lenVertexProc;
      float deltaMax = N*(epsilon -1);

      int modProc = N - numProcs*lenVertex;

      if (procRank < modProc){
        lenVertexProc = lenVertex+1;
      }else{
        lenVertexProc = lenVertex;
      }

      int *nextPart = malloc(sizeof(part[0])*lenVertexProc);

      // array used for allgatherv
      int *disp = malloc(sizeof(int)*numProcs);
      int *recvCount = malloc(sizeof(int)*numProcs);
      for (int i = 0; i<numProcs; i++){
        if (i<modProc){
          disp[i]= i * (1+lenVertex);
          recvCount[i]=1+lenVertex;
        }else{
          disp[i]= i * lenVertex+modProc;
          recvCount[i] = lenVertex;
        }
      }


      deltaTot = calculDeltaPart(N,part);

      for (int i = 1; i<= maxIters; i++){
        for (int j = 0; j< lenVertexProc; j++){
          currentVertex = offsetProc + j;
          signPart = convert(part[currentVertex]);
          if (-((deltaTot+deltaLocal) * signPart) > deltaMax){
            nextPart[j] = part[currentVertex];
          } else{
            nextPart[j] = futurLabel(currentVertex, adjBeg,adj, part);
            if (nextPart[j] != part[currentVertex]){
              deltaLocal -= 2*signPart;
            }
          }
        }

        // communications

        MPI_Allgatherv(nextPart, lenVertexProc, MPI_INT,
                       part, recvCount, disp, MPI_INT, MPI_COMM_WORLD);
        MPI_Allreduce(&deltaLocal, &nextDelta, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        deltaTot += nextDelta;
        nextDelta = 0;
        deltaLocal = 0;

        printCutSizeImbal(N, adjBeg, adj, part);
      }
      free(nextPart);
      free(disp);
      free(recvCount);



  } else if (strcmp(algoName, "bisect-p2p") == 0) { // Graph bisection with point-to-point communication
    //initialization

    int signPart;
    int deltaTot = 0;
    int nextDelta = 0;
    int deltaLocal = 0;
    float deltaMax = N*(epsilon -1);
    int totalVertexSend;
    int totalVertexRecv;

    //general variables similar to a2a
    int lenVertex = N /numProcs;
    int offsetProc = lenVertex*procRank;
    int currentVertex;
    int lenVertexProc;
    int modProc = N - numProcs*lenVertex;

    if (procRank < modProc){
      lenVertexProc = lenVertex+1;
      offsetProc= procRank * (1+lenVertex);
    }else{
      lenVertexProc = lenVertex;
      offsetProc= procRank * lenVertex+modProc;
    }

    int *nextPart = malloc(sizeof(int)*N);

    // initialisation : computation of the number of elements to communicate
    int *recvCount = calloc(sizeof(int),numProcs);
    int *vertexNeeded = calloc(sizeof(int),N);

    for (int v = offsetProc; v<lenVertexProc+offsetProc;v++){
      detectNeighbour(N, adjBeg, adj, v, vertexNeeded, recvCount, numProcs);
    }

    int *recvIdxBeg = malloc(sizeof(int)*(numProcs+1));
    cumulativeSum(recvCount, recvIdxBeg, numProcs);
    totalVertexRecv = recvIdxBeg[numProcs];

    // creation of the array with the computed sizes
    memset(vertexNeeded, 0, sizeof(int)*N);
    int *recvIdx = malloc(sizeof(int)*totalVertexRecv);
    int *currentIdx = calloc(sizeof(int),numProcs); // aux array
    for (int v = offsetProc; v<lenVertexProc+offsetProc;v++){
      assignNeighbour(N, adjBeg, adj, v, vertexNeeded, currentIdx, recvIdx, recvIdxBeg, numProcs);
    }
    // initialization for the sending

    // similarly to send, compute the size
    int *sendCount = calloc(sizeof(int),numProcs);
    // communicate the size
    MPI_Alltoall(recvCount, 1, MPI_INT, sendCount, 1, MPI_INT, MPI_COMM_WORLD);

    //create arrays thanks to the given information
    int *sendIdxBeg = malloc(sizeof(int)*(numProcs+1));
    cumulativeSum(sendCount, sendIdxBeg, numProcs);
    totalVertexSend = sendIdxBeg[numProcs];
    int *sendIdx = malloc(sizeof(int)*totalVertexSend);

    //compute the array
    MPI_Alltoallv(recvIdx, recvCount, recvIdxBeg, MPI_INT, sendIdx, sendCount, sendIdxBeg, MPI_INT, MPI_COMM_WORLD);

    // initialization of the delta
    deltaTot = calculDeltaPart(N,part);

    int *sendPart = malloc(sizeof(int)*totalVertexSend);
    int *recvPart = malloc(sizeof(int)*recvIdxBeg[numProcs]);

    for (int i = 1; i<= maxIters; i++){
      for (int j = 0; j< lenVertexProc; j++){
        currentVertex = offsetProc + j;
        signPart = convert(part[currentVertex]);
        if (-((deltaTot+deltaLocal) * signPart) > deltaMax){
          nextPart[currentVertex] = part[currentVertex];
        } else{
          nextPart[currentVertex] = futurLabel(currentVertex, adjBeg,adj, part);
          if (nextPart[currentVertex] != part[currentVertex]){
            deltaLocal -= 2*signPart;
          }
        }
      }

      // communications

      for(int v = 0 ; v < totalVertexSend; v++){
        sendPart[v] = nextPart[sendIdx[v]];
      }

      MPI_Alltoallv(sendPart, sendCount, sendIdxBeg, MPI_INT, recvPart, recvCount, recvIdxBeg, MPI_INT, MPI_COMM_WORLD);


      for (int v = 0 ; v<totalVertexRecv; v++){
        part[recvIdx[v]] = recvPart[v];
      }

      // same as a2a : pass on the delta
      MPI_Allreduce(&deltaLocal, &nextDelta, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      deltaTot += nextDelta;
      nextDelta = 0;
      deltaLocal = 0;


      printCutSizeImbal(N, adjBeg, adj, part);
    }

    free(recvIdx);
    free(sendIdx);
    free(recvIdxBeg);
    free(sendIdxBeg);
    free(vertexNeeded);
    free(currentIdx);
    free(recvCount);
    free(sendCount);
    free(nextPart);
    free(sendPart);
    free(recvPart);

  } else { // Invalid algorithm name
    if (procRank == 0) {
      printf("ERROR: Invalid algorithm name: %s.\n", algoName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}