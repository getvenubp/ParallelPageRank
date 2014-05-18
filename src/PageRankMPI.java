import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import mpi.*;

/**
 * @author Venu
 */
public class PageRankMPI {

	static HashMap<Integer, ArrayList<Integer>> adjacencyMatrix = new HashMap<Integer, ArrayList<Integer>>();
	static int[] amIndex=null;
	static int[] adjMatrix=null;

	/**
	 * MPI Read will read from the file and will split it among processors
	 * @param inputFileName
	 * @param numUrls
	 * @param comm
	 * @return localNumUrls
	 */
	public static int mpiRead(String inputFileName,int numUrls,Intracomm comm)throws IOException{
		
		//getting rank of the processor
		int rank =  comm.Rank();
		//getting total number of processors
		int nproc = comm.Size();

		int localNumUrls=0;
		int totalNumUrls=0;
		int divisions = 0;
		int remainder = 0;
		int sendStartIndex = 0;

		int[] totalUrls= new int[1];
		//reading from input file specified and creating Adjacency Matrix
		if(rank==0)
			totalUrls[0]=fileRead(inputFileName,totalNumUrls);
		
		//broadcasting total number of urls
		comm.Bcast(totalUrls, 0, 1, MPI.INT, 0);		
		totalNumUrls=totalUrls[0];
		
		if(rank==0){
			
			//dividing total number of urls by processors to determine blocks
			divisions = totalNumUrls / nproc;
			remainder = totalNumUrls % nproc;

			int k=0;
			
			for(int i=0;i<nproc;i++){
				sendStartIndex=k;
				
				if(remainder == 0)
					numUrls = divisions;
				else
					numUrls = (i < remainder) ? (divisions + 1) : divisions;
				
				int index=0;
				int len = 0;
				int[] sendamIndex = new int[numUrls*2];
				
				//creating am index array to send to processors
				for(k = sendStartIndex; k <= ((sendStartIndex+numUrls)-1) ; k++){
					sendamIndex[index++] = k;
					sendamIndex[index++] = (adjacencyMatrix.get(k).size() + 1);
					len+=adjacencyMatrix.get(k).size()+1;
				}

				index=0;
				int[] sendadjMatrix = new int[len];
				//creating adjMatrix array to send to processors
				for(k = sendStartIndex; k <= ((sendStartIndex+numUrls)-1) ; k++){
					sendadjMatrix[index++] = k;
					Iterator<Integer> value = adjacencyMatrix.get(k).iterator();
					while(value.hasNext()){
						sendadjMatrix[index++] = value.next();
					}
				}
				 
				if(i==0){
					//assigning am index array, adjacencyMatrix array and number of urls to rank 0
					amIndex = new int[numUrls*2];
					adjMatrix = new int[len];
					amIndex = sendamIndex;
					adjMatrix = sendadjMatrix;
					localNumUrls=numUrls;
				}
				else{
					int[] sendamIndexSize = new int[1];
					sendamIndexSize[0] = numUrls*2;
					 
					//sending am index array size
					comm.Send(sendamIndexSize, 0, 1, MPI.INT, i, 0);
					//sending am index array
					comm.Send(sendamIndex, 0, sendamIndexSize[0], MPI.INT, i, 1);
	                 
	                int[] sendadjMatrixSize = new int[1];
	                sendadjMatrixSize[0] = len;
	          	   	
	                //sending adjacencyMatrix array size
	                comm.Send(sendadjMatrixSize, 0, 1, MPI.INT, i, 2);
	              	//sending adjacencyMatrix array
	                comm.Send(sendadjMatrix, 0, sendadjMatrixSize[0], MPI.INT, i, 3);
	            }
			}
		}// end of if rank==0
		else{
			int[] recvamIndexSize= new int[1];
			//Receiving am index array size
	        comm.Recv(recvamIndexSize, 0, 1, MPI.INT, 0, 0);
	        
	        localNumUrls= recvamIndexSize[0]/2;
	        
	        amIndex = new int[recvamIndexSize[0]];
	        //Receiving am index array
	        comm.Recv(amIndex, 0, recvamIndexSize[0], MPI.INT, 0, 1);

	        int[] adjMatrixSize = new int[1];
	      	//Receiving adjacencyMatrix array size
	        comm.Recv(adjMatrixSize, 0, 1, MPI.INT, 0, 2);
	        
	        adjMatrix = new int[adjMatrixSize[0]];
	      	//Receiving adjacencyMatrix array
	        comm.Recv(adjMatrix, 0, adjMatrixSize[0], MPI.INT, 0, 3);
		
		}// end of else
	return localNumUrls;
	}// end of mpiRead function
	
	/**
	 * Parses the input file and initializes the adjacency matrix
	 * @param inputFileName
	 * @param totalNumUrls
	 * @throws IOException
	 */
	public static int fileRead(String inputFileName,int totalNumUrls) throws IOException{
	
		FileReader fr=null;
		BufferedReader br = null;
		try{
			//Opening the file specified to read
			fr= new FileReader(inputFileName);
			br = new BufferedReader(fr);
			
			String str;
			Integer nodeNum = null;
			
			//Reading one line at a time from file and adding it hash map
			while((str=br.readLine())!=null){
				ArrayList <Integer> intermediateList = new ArrayList <Integer>();
				String stringRead[]=str.split(" ");
				nodeNum = Integer.parseInt(stringRead[0]);
				for(int i=1;i<stringRead.length;i++){
					intermediateList.add(Integer.parseInt(stringRead[i]));
				}
				adjacencyMatrix.put(nodeNum,intermediateList);
			}
			totalNumUrls=adjacencyMatrix.size();
		}
		//catching file not found exception
		catch(FileNotFoundException fe){
			System.out.println("File Not Found Exception");
			fe.printStackTrace();
		}
		//catching any other exception
		catch(Exception e){
			System.out.println("Exception in opening a file");
			e.printStackTrace();
		}
		//at the end closing file reader and buffer reader
		finally{
			br.close();
			fr.close();
		}
	return totalNumUrls;
	}// end of fileRead

	/**
	 * Assigns initial rank values for each page in rankValuesTable as 1/N
	 * @param rankValuesTable
	 * @param totalNumUrls
	 */
	public static void initialRankTable(double[] rankValuesTable, int totalNumUrls){
	   
	   double initialRankValuePerPage = 1.0/(double)totalNumUrls;
	   
	 //assigns intital values to rankValuesTable
	   for (int i=0;i<totalNumUrls;i++)
			rankValuesTable[i] = initialRankValuePerPage;
	}

	/**
	 * Calculates mpiPageRank
	 * @param numUrls
	 * @param totalNumUrls
	 * @param iterations
	 * @param threshold
	 * @param rankValuesTable
	 * @param comm
	 */
	public static void mpiPageRank(int numUrls,int totalNumUrls, int iterations, double threshold, double[] rankValuesTable,Intracomm comm){  
		
		double[] localRankValuesTable = new double[totalNumUrls];
		double[] oldRankValuesTable = new double[totalNumUrls];
		double[] delta = new double[1];
		int i, j, k=1;
		//getting rank of the processor
		int rank = comm.Rank();

		//assigning initial value of rankvalueTable to oldRankValuesTable
		for(i=0;i<totalNumUrls;i++)
			oldRankValuesTable[i] = rankValuesTable[i];

	do{
		int sourceURL=0,targetURL=0,index=0;
		double value =0.0,dangling = 0.0;
		double dampingFactor = 0.85;
		
		for(i=0;i<totalNumUrls;i++)
			localRankValuesTable[i] = 0.0;
		
		for(i=1;i<amIndex.length;i=i+2){
			sourceURL=adjMatrix[index++];
			value =rankValuesTable[sourceURL]/(double)(amIndex[i]-1);
			for(j=1;j<amIndex[i];j++){
				targetURL=adjMatrix[index++];
				localRankValuesTable[targetURL]+=value;			
			}
			if(amIndex[i]==1)
				dangling += rankValuesTable[sourceURL];
		}
		
		//getting localRankValuesTable from all processors
		comm.Allreduce(localRankValuesTable, 0, rankValuesTable, 0, totalNumUrls, MPI.DOUBLE, MPI.SUM);
					
		double[] danglingValue = new double[1];
		danglingValue[0] = dangling;
		double[] danglingValueSum = new double[1];
		
		comm.Allreduce(danglingValue, 0, danglingValueSum, 0, 1, MPI.DOUBLE, MPI.SUM);
		
		if(rank==0){
			double sum=0.0;
			double avgDanglingValue = danglingValueSum[0]/totalNumUrls;
		
			//Computing page rank
			for(i=0;i<totalNumUrls;i++){
				rankValuesTable[i] = rankValuesTable[i] + avgDanglingValue;
				rankValuesTable[i] = dampingFactor * rankValuesTable[i] + (1 - dampingFactor) * (1.0/(double)totalNumUrls);
			}
			
			for(i=0;i<totalNumUrls;i++)
				sum=sum+rankValuesTable[i];
			
			delta[0] = 0.0;
			double diff = 0.0;
			for(i=0;i<totalNumUrls;i++){
				diff = oldRankValuesTable[i] - rankValuesTable[i];
				delta[0] += diff * diff; 
				oldRankValuesTable[i] = rankValuesTable[i];
			}
			System.out.println("Current iteration = "+k+" Delta = "+delta[0]);
		}
		//broadcoasting delta to other processors
		comm.Bcast(delta, 0, 1, MPI.DOUBLE, 0);
		//broadcasting rankValueTable to other processors
		comm.Bcast(rankValuesTable, 0, totalNumUrls, MPI.DOUBLE, 0);
		
	}while(k++<iterations&&delta[0]>threshold);
}

	/**
	 * Writes top 10 urls into an output file by Sorts urls according to page rank value  
	 * @param rankValueTable
	 * @param numOfUrls
	 * @throws IOException
	 */
	public static void mpiWrite(String outputFileName,int totalNumUrls,double[] rankValuesTable) throws IOException{
		LinkedHashMap <Integer,Double> sortedPageRank = new LinkedHashMap <Integer,Double>();
		HashMap <Integer,Double> rankTable = new HashMap <Integer,Double>();
		
		for(int i = 0;i < totalNumUrls;i++)
			rankTable.put(i,rankValuesTable[i]);
		
		ArrayList<Integer> rankValuesTableKeys = new ArrayList<Integer>(rankTable.keySet());
		ArrayList<Double> rankValuesTableValues = new ArrayList<Double>(rankTable.values());
		
		//sorting all rank value is descending order
		Collections.sort(rankValuesTableValues,Collections.reverseOrder());
		Iterator<Double> valueIterator = rankValuesTableValues.iterator();
		
		//adding sorted rank values to linked hash map
		while(valueIterator.hasNext()){
			Double value = valueIterator.next();
			Iterator<Integer> keyIterator = rankValuesTableKeys.iterator();
			while(keyIterator.hasNext()){
				Integer key = keyIterator.next();
				Double pageRankValue = rankTable.get(key);
				if(value.equals(pageRankValue)){
					rankTable.remove(key);
					rankValuesTableKeys.remove(key);
					sortedPageRank.put(key, value);
					break;
				}
			}
		}
		
		
		FileWriter file=null;
		BufferedWriter out =null;
		try{
			file = new FileWriter(outputFileName);
			out = new BufferedWriter(file);
			ArrayList<Integer> sortedPageRankKeys = new ArrayList<Integer>(sortedPageRank.keySet());
			Iterator<Integer> iterator = sortedPageRankKeys.iterator();
			out.write("Top 10 Page Ranked URLs are:\n");
			out.write("----------------------------\n");
			out.write("URL\tPageRank value\n");
			out.write("----------------------------\n");
			int i=0;
			double sum=0.0;
			
			//writing top 10 URLs to file
			while(iterator.hasNext()){
				Integer key = iterator.next();
				if(i<10){
					out.write(key+"\t"+sortedPageRank.get(key)+"\n");
					i++;
				}
				sum+=sortedPageRank.get(key);
			}
			out.write("----------------------------\n");
			System.out.println("Sum of all page rank = "+sum);
		}
		//catching file not found exception
		catch(FileNotFoundException fe){
			System.out.println("File Not Found Exception\n");
			fe.printStackTrace();
		}
		//catching other exceptions
		catch(Exception e){
			System.out.println("Exception caught\n");
			e.printStackTrace();
		}
		//closing file write and buffer writer
		finally{
			out.close();
			file.close();
		}
		
	}// end of mpiWrite function

	/**
	 * Main Method 
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		
		if(args.length!=7){
		System.out.println("\njava PageRankMPI [input filename] [output filename] [iteration count] [threshold]");
		System.exit(-1);
		}
		
		String inputFileName = args[3];
		String outputFileName = args[4];
		int iterations = Integer.parseInt(args[5]);
		double threshold = Double.parseDouble(args[6]);
		
		if(iterations<0){
		System.out.println("\n[iteration count] cannot be less than zero");
		System.exit(-1);
		}
		
		if(threshold<=0 || threshold>=1){
		System.out.println("\n[threshold] should be between 0 and 1");
		System.exit(-1);
		}

		int numUrls=0;
		int totalNumUrls=0;
		int[] numOfUrls= new int[1];
		int[] totalUrls= new int[1];

		MPI.Init(args);

		long start1 = System.currentTimeMillis();
		
		int rank =  MPI.COMM_WORLD.Rank();
		int nproc = MPI.COMM_WORLD.Size();
		Intracomm comm = MPI.COMM_WORLD;
		
		//calling mpiRead
		numOfUrls[0]=mpiRead(inputFileName,numUrls,comm);
		numUrls=numOfUrls[0];
		
		//getting total number of urls 
		comm.Allreduce(numOfUrls, 0, totalUrls, 0, 1, MPI.INT, MPI.SUM);
		totalNumUrls=totalUrls[0];
		
		double[] rankValuesTable= new double[totalNumUrls];

		if(rank==0)
			//calling initialRankTable to initialize rank value table
			initialRankTable(rankValuesTable,totalNumUrls);
		
		//broadcasting rank value table to other processors
		comm.Bcast(rankValuesTable, 0, totalNumUrls, MPI.DOUBLE, 0);		

		if(rank==0)
			System.out.println("Max iterations = "+iterations+" Threshold = "+threshold);

		long start2 = System.currentTimeMillis();
		
		//calling mpiPageRank to compute page rank
		mpiPageRank(numUrls,totalNumUrls,iterations,threshold,rankValuesTable,comm);
		
		long finish2 = System.currentTimeMillis();

		if(rank==0){
			//calling mpiWrite to sort and write top 10 urls to file
			mpiWrite(outputFileName,totalNumUrls,rankValuesTable);
			long finish1 = System.currentTimeMillis();
			System.out.println("Proc:"+rank+" has written rank values of 10 urls to file "+outputFileName);
			System.out.println("*****MPI Page Rank******");
			System.out.println("Number of processes = "+nproc);
			System.out.println("Input filename = "+inputFileName);
			System.out.println("Total number of urls = "+totalNumUrls);
			System.out.println("Number of iterations = "+iterations);
			System.out.println("Threshold = "+threshold);
			System.out.println("I/O time = "+(finish1 - start1)/1000.0+" sec");
			System.out.println("Computation timing = "+(finish2 - start2)/1000.0+" sec");
		}
		
		MPI.Finalize();

	}// end of main function

}//end of class PageRankMPI
