/***************************************************************************

 Object version :
    Sang Lim(slim@npac.syr.edu)
    Northeast Parallel Architectures Center at Syracuse University
    11/15/98
****************************************************************************/

import mpi.*;
 
class gatherO {
  static public void main(String[] args) throws MPIException {
    
    int root=0,i,j,k,l;
    int myself,tasks;
 
    MPI.Init(args);
    myself = MPI.COMM_WORLD.Rank();
    tasks = MPI.COMM_WORLD.Size();
 
    int out[][] = new int[6][3];
    int in[][]  = new int[6*tasks][];

    for (j = 0; j < tasks; j++){
      if (j == myself)
        for (l=0;l<6;l++)
          for (i = 0; i < 3; i++)
            out[l][i] = l+j*6;
    }

    
    MPI.COMM_WORLD.Gather(out,0,6,MPI.OBJECT,in,0,6,MPI.OBJECT,root);

    if (root == myself){
      for (l=0;l<6*tasks;l++)
	for (i = 0; i < 3; i++)
          if (in[l][i] != l)
            System.out.println("Recived data : "+in[l][i]+" at recive buffer["+l+"]["+i+"] should be : "+ l);
    }

    MPI.COMM_WORLD.Barrier();
    if(myself == root)  System.out.println("GatherO TEST COMPLETE\n");
    MPI.Finalize();
  }
}
