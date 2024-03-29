import mpi.*;

class bsend {
  static public void main(String[] args) throws MPIException {
    /* Note that the buffer sizes must include the BSEND_OVERHEAD;
       these values are probably sizeof(int) too large */
    int len,tsks,me,i,size,rc;	
    final int A100 = 1000000;
    Status status;

    MPI.Init(args);
    me=MPI.COMM_WORLD.Rank();

    int data100[] = new int[A100];
    //int intsize = 4;
    int intsize = MPI.INT.Size();

    byte buf100[] = new byte[A100*intsize+MPI.BSEND_OVERHEAD];

    if ( me == 0 ) {      
      MPI.Buffer_attach(buf100);
      /* test to see if large array is REALLY being buffered */
      for(i=0;i<A100;i++)  data100[i] = 1;
      
      MPI.COMM_WORLD.Bsend(data100,0,A100,MPI.INT,1,1);
      
      MPI.COMM_WORLD.Recv(data100,0,A100,MPI.INT,1,2);
      
      for(i=0;i<A100;i++)
	if(data100[i] != 2)  
	  System.out.println
	    ("ERROR, incorrect data["+i+"]="+data100[i]+", task 0");
    } else if ( me == 1 ) {
      MPI.Buffer_attach(buf100);
      /* test to see if large array is REALLY being buffered */
      for(i=0;i<A100;i++)  data100[i] = 2;
      
      MPI.COMM_WORLD.Bsend(data100,0,A100,MPI.INT,0,2);      
      MPI.COMM_WORLD.Recv(data100,0,A100,MPI.INT,0,1);

      for(i=0;i<A100;i++)
	if(data100[i] != 1)  
	  System.out.println
	    ("ERROR , incorrect data["+i+"]="+data100[i]+", task 1");
      
    }

    MPI.COMM_WORLD.Barrier();
    if(me == 0)  System.out.println("Bsend TEST COMPLETE\n");
    MPI.Finalize();
  }

}
