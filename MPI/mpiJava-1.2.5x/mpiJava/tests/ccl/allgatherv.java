/*
 MPI-Java version :
    Sang Lim (slim@npac.syr.edu)
    Northeast Parallel Architectures Center at Syracuse University
    12/1/98
*/

import mpi.*;
 
class allgatherv {
  static public void main(String[] args) throws MPIException {
    
    final int MAXLEN = 10;
 
    int root,i,j,k;
    int myself,tasks,stride=15;
 
    MPI.Init(args);
    myself = MPI.COMM_WORLD.Rank();
    tasks = MPI.COMM_WORLD.Size();

    int out[] = new int[MAXLEN];
    int in[]  = new int[MAXLEN*stride*tasks];
    int dis[] = new int[MAXLEN];
    int rcount[] = new int[MAXLEN];
    int ans[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    for (i = 0; i<MAXLEN;i++){
       dis[i] = i*stride;
       rcount[i] = 5;
       out[i] = i;
    }
    rcount[0] = 10;

    for (i = 0; i<MAXLEN*tasks*stride;i++){
       in[i] = 0;
    }

    if (myself == 0)
       MPI.COMM_WORLD.Allgatherv(out,0,10,MPI.INT,
                              in ,0,rcount,dis,MPI.INT);
    else 
       MPI.COMM_WORLD.Allgatherv(out,0,5,MPI.INT,
                              in ,0,rcount,dis,MPI.INT);
    
    for(i=0; i<tasks*stride; i++)
      if (ans[i]!=in[i])
        System.out.println("recived data : "+in[i]+"at ["+i+"] should be : "+ans[i]+" on proc. : "+myself);

      MPI.COMM_WORLD.Barrier();

    if(myself == 0)  System.out.println("Allgatherv TEST COMPLETE\n");
    MPI.Finalize();
  }
}
