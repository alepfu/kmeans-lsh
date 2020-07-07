/*
 * Used references:
 *	Time measurement
 *		http://stackoverflow.com/questions/1770010/how-do-i-measure-time-elapsed-in-java
 *	CSV file parsing
 *		http://stackoverflow.com/questions/6857248/fast-csv-parsing
 *		http://javarevisited.blogspot.co.at/2016/01/reading-writing-files-using-filechannel-bytebuffer-example.html
 *  Saving a buffer as CSV file
 *      http://stackoverflow.com/questions/735600/how-to-save-a-java-floatbuffer-or-any-buffer-to-a-file
 *	Manual memory management: java.nio -> FloatBuffer
 *		http://www.javamex.com/java_equivalents/malloc.shtml
 *		http://stackoverflow.com/questions/10697161/why-floatbuffer-instead-of-float
 *		https://docs.oracle.com/javase/8/docs/api/java/nio/FloatBuffer.html
 *		https://www.ntu.edu.sg/home/ehchua/programming/java/J5b_IO_advanced.html
 *      http://stackoverflow.com/questions/28744096/convert-bytebuffer-to-byte-array-java
 *	sqrt() of float numbers
 *  	http://stackoverflow.com/questions/15143735/using-floatmath-or-math-and-a-cast
 */

package at.alepfu;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;


public class Main {

	/*
	   Input parameters
	   ----------------
	*/
	final static String filename = "LSH data.txt";
	final static String outputFilename = "kmeans_result.txt";
	final static String bucketSizeDistributionFilename = "bucket_distribution.txt";
	final static int NUM_POINTS = 291500;
	final static int NUM_DIM = 10;
	final static int NUM_CLUSTERS = 15;

	final static int NUM_HASHFUNCTIONS = 25;    // max. 30
	final static int NUM_BUCKETS = 33554432;	// must equal 2^NUM_HASHFUNCTIONS

	// HASH_SCALING influences the proportion of filled buckets
	// (increased scaling factor -> more filled buckets)
	final static float HASH_SCALING = 1f;

	// other operational parameters
	final static int MAX_ITER = 500;	// maximum number of k-means iterations
	final static int SEED = 137;		// for reproducibility

	/*
	   DO NOT EDIT BELOW
	   -----------------
	*/
	// operational helpers
	final static int EMPTY_BUCKET = -1;

	// RNG
	final static Random rand = new Random(SEED);
	//final static Random rand = new Random();

	// hash functions
	// NB: type float is used throughout because it gives a few percent
	//     improvement in runtime performance (vs. double)
	static float[][] randvect = new float[NUM_HASHFUNCTIONS][NUM_DIM];
	static float[] scaling = new float[NUM_HASHFUNCTIONS];

	/*
		Buckets
		-------
		* Each point belongs to exactly one bucket.
		* A bucket is defined as the combination of all (binary) hash functions
		  and is represented as 32bit integer.
	*/
	// for each point: what is its bucket
	final static int[] bucketIdxOfPoint = new int[NUM_POINTS];
	// map of bucket to representative (index):
	//   for each bucket: what is its representative point
	//   (index into 1..NUM_POINTS)
	final static HashMap<Integer,Integer> bucket2RepPointIdxMap = new HashMap<>();
	// map from a bucket (hash value) to a linear index (1..# of filled buckets)
	final static HashMap<Integer,Integer> bucket2bucketIdxMap = new HashMap<>();


	/**
	 * Initialize hash functions randomly.
	 * Random Gaussians are created by adding 12 uniformly distributed integers.
	 * We expect these to perform sufficiently in our application and also much
	 * faster than java.util.nextGaussian, the latter using Box-Muller algo.
	 * (https://docs.oracle.com/javase/7/docs/api/java/util/Random.html#nextGaussian%28%29)
	 */
	public static void init_hashfunctions() {
		// random vector
		float sumrand;
		for( int h = 0; h < NUM_HASHFUNCTIONS; h++ ) {
			// generate random number vector
			for (int m = 0; m < NUM_DIM; m++) {
				sumrand = 0.0f;
				for (int x = 0; x < 12; x++) {
					sumrand += rand.nextFloat();
				}
				randvect[h][m] = sumrand - 6;
			}
			// scaling parameters (similar to w parameter on lecture slides)
			scaling[h] = HASH_SCALING;
		}
	}


	/**
	 * Initialize (hash) buckets from (original space) points.
	 * @param points
	 * @return number of non-empty buckets
	 */
	public static int points2buckets( final float[][] points ) {

		float dotprod;
		int bucket;
		int bucketIdx = 0;

		/*
			put points into buckets & select first point of each bucket as its
			representative
		*/
		for (int n = 0; n < NUM_POINTS; n++) {
			bucket = 0;
			for( int h = 0; h < NUM_HASHFUNCTIONS; h++ ) {
				dotprod = 0.0f;
				for (int m = 0; m < NUM_DIM; m++) {
					dotprod += points[n][m] * randvect[h][m];
				}
				dotprod *= scaling[h];

				// way of bucketing with sign (results in suboptimal bucketing)
				//final int sign = ( dotprod < 0f ? 1 : 0 );

				// bucketing with sign & modulo 256
				// http://stackoverflow.com/a/20369990
				bucket = ( bucket << 1 ) + (new Float(dotprod).byteValue() >>> 31 );
			}
			assert bucket >= 0;
			assert bucket <= Integer.MAX_VALUE;

			// select first point as representative
			final Integer existingBucketIdx = bucket2bucketIdxMap.get(bucket);
			if( null == existingBucketIdx ) {
				assert !bucket2RepPointIdxMap.containsKey(bucket);
				bucket2RepPointIdxMap.put(bucket, n);
				bucket2bucketIdxMap.put(bucket, bucketIdx);
				// conserve mapping of point to bucket index
				bucketIdxOfPoint[n] = bucketIdx++;
			} else {
				// conserve mapping of point to bucket index
				bucketIdxOfPoint[n] = existingBucketIdx;
			}

		}

		// return number of filled (i.e. non-empty) buckets
		return bucketIdx;
	}


	/**
	 * Read points from file.
	 * OUTPUTS:
	 * @param points
	 * @throws IOException
	 */
	public static void read_points_from_file( final float[][] points ) throws IOException {

		String line;
		BufferedReader reader;
		int i, j;

		reader = new BufferedReader(new FileReader(filename)); // buffer size can be specified via extra argument
		i = 0;
		while ((line = reader.readLine()) != null) {
			j = 0;
			for(final String entry : line.split(",")) {
				points[i][j++] = Float.parseFloat(entry);
			}
			i++;
		}
		reader.close();
	}


	/**
	 * Save results together with input points (in original data format) to file.
	 * @param pointBelongsToCluster: the clustering result
	 * @throws IOException
	 */
	public static void write_results_to_file( final int[] pointBelongsToCluster ) throws IOException {

		final BufferedReader reader = new BufferedReader(new FileReader(filename));
		final FileWriter writer = new FileWriter(outputFilename);

		int lineIdx = 0;
		String line;

		while ((line = reader.readLine()) != null) {
			writer.append( line );
			writer.append( ',' );
			writer.append( Integer.toString(pointBelongsToCluster[lineIdx++]) );
			writer.append( System.lineSeparator() );
		}
		reader.close();
		writer.close();
	}


	/**
	 * Print analysis of the clustering given by
	 * @param points
	 * @param pointBelongsToCluster
	 * to STDOUT, while freshly calculating all relevant statistics from scratch.
	 */
	public static void print_clustering_analysis( final float[][] points, final int[] pointBelongsToCluster, final String type ) {
		final int NUM_POINTS = points.length;
		final int NUM_DIM = points[0].length;
		int NUM_CLUSTERS;
		String cluster_string;

		if( type == "bucket" ) {
			NUM_CLUSTERS = bucket2bucketIdxMap.size();
			cluster_string = "Bucket";
		} else {
			NUM_CLUSTERS = Main.NUM_CLUSTERS;
			cluster_string = "Cluster";
		}

		final float[][] centroids = new float[NUM_CLUSTERS][NUM_DIM];
		final double[] gyration_radii = new double[NUM_CLUSTERS];
		double mean_radius = 0;
		final int[] cluster_size = new int[NUM_CLUSTERS];
		float compactness = 0;

		// calculate cluster sizes, centroids
		for( int p = 0; p < NUM_POINTS; p++ ) {
			final int k = pointBelongsToCluster[p];
			for( int d = 0; d < NUM_DIM; d++ )
				centroids[k][d] += points[p][d];
			cluster_size[k]++;
		}
		for( int k = 0; k < NUM_CLUSTERS; k++ )
			for( int d = 0; d < NUM_DIM; d++ )
				centroids[k][d] /= cluster_size[k];

		// calculate gyration radii & compactness
		for( int p = 0; p < NUM_POINTS; p++ ) {
			final int k = pointBelongsToCluster[p];
			for( int d = 0; d < NUM_DIM; d++ ) {
				final float dist = centroids[k][d] - points[p][d];
				gyration_radii[k] += dist*dist;
			}
		}
		for( int k = 0; k < NUM_CLUSTERS; k++ ) {
			compactness += gyration_radii[k];
			if( gyration_radii[k] > 0 ) {
				gyration_radii[k] = Math.sqrt( gyration_radii[k] / cluster_size[k] );
				mean_radius += gyration_radii[k] * cluster_size[k];
			}
		}
		mean_radius /= NUM_POINTS;

		// print results
		System.out.println(String.format("%sing compactness: TD2 = %.6g [per point %.6g]", cluster_string, compactness, compactness/NUM_POINTS ));
		System.out.println(String.format("Weighted mean gyration radius <R> = %.6g", mean_radius ));

		if( type != "bucket" ) {
			for( int k = 0; k < NUM_CLUSTERS; k++ ) {
				if( 0 == cluster_size[k] && type=="bucket" ) continue;
				System.out.print(String.format("%s %2d [%6d points, R=%4.1f]: ", cluster_string, k, cluster_size[k], gyration_radii[k]));
				for( int d = 0; d < NUM_DIM; d++ )
					System.out.print(String.format("%6.2f ", centroids[k][d]));
				System.out.println();
			}
		}

		// save analysis to file
		if( type == "bucket" ) {
			FileWriter fileWriter;
			try {
				fileWriter = new FileWriter(bucketSizeDistributionFilename);
				for( int k = 0; k < NUM_CLUSTERS; k++ ) {
					fileWriter.append( new Integer(cluster_size[k]).toString() );
					fileWriter.append( "\t" );
					fileWriter.append( new Double(gyration_radii[k]).toString() );
					fileWriter.append( System.lineSeparator() );
				}
				fileWriter.close();
			} catch (final IOException e) {
				e.printStackTrace();
			}
		}
	}


	/**
	 * Perform k-means clustering on points
	 * using I-1 (random points as centroids) and U-2 (MacQueen).
	 * INPUTS:
	 * @param points : input points
	 * OUTPUTS:
	 * @param centroids: centroid coordinates
	 * @param clusterContainsNumPoints: number of points per cluster
	 * @param pointBelongsToCluster: cluster memberships
	 * @return number of iterations
	 */
	public static int kmeans(
			final float[][] points,
			final float[][] centroids,
			final int[] clusterContainsNumPoints,
			final int[] pointBelongsToCluster
			) {

		final int NUM_POINTS = points.length;
		final int NUM_DIM = points[0].length;
		final int NUM_CLUSTERS = centroids.length;

		// validate inputs (needs VM argument '-ea')
		assert NUM_DIM == centroids[0].length;
		assert NUM_CLUSTERS == clusterContainsNumPoints.length;
		assert NUM_POINTS == pointBelongsToCluster.length;

		// Declare all variables
		int iter, i, j, k;
		int previouslyAssignedCluster, newlyAssignedCluster;
		boolean membershipChange;
		float diff, currentSqDist, smallestSqDist;
		int nearestCentroid = 0;
		int N_previouslyAssignedCluster;
		int N_newlyAssignedCluster;

		// declare & allocate & initialize (zeros by Java specs) arrays
		final int [] randomEntries = new int[NUM_CLUSTERS];


		/*
			1) Initialize clusters.
			Initialization strategy I-1: Random centroids.
		*/

		// randomly select points to be used as initial cluster centers
		for(i=0; i<NUM_CLUSTERS; i++)
			randomEntries[i] = rand.nextInt(NUM_POINTS); // TODO: guarantee uniqueness of initial cluster centers

		// Choose random centroids as cluster centers
		for(i=0; i<NUM_CLUSTERS; i++)
			for(j=0; j<NUM_DIM; j++)
				centroids[i][j] = points[ randomEntries[i] ][j];

		// Assign each point to the cluster whose center is in shortest distance from the point
		for(i=0; i<NUM_POINTS; i++) {
			smallestSqDist = Float.MAX_VALUE;
			for(k=0; k<NUM_CLUSTERS; k++) {
				currentSqDist = 0.0f;
				for(j=0; j<NUM_DIM; j++) {
					diff = points[i][j] - centroids[k][j];
					currentSqDist += diff * diff;
				}

				if(currentSqDist < smallestSqDist) {
					smallestSqDist = currentSqDist;
					nearestCentroid = k;
				}
			}
			pointBelongsToCluster[i] = nearestCentroid;
			clusterContainsNumPoints[nearestCentroid]++;
		}

		// Re-determine all cluster centers
		// Reset cluster centers to 0.0
		for(k=0; k<NUM_CLUSTERS; k++)
			for(j=0; j<NUM_DIM; j++)
				centroids[k][j] = 0.0f;

		// Calculate the sum over all points in each cluster
		for(i=0; i<NUM_POINTS; i++)
			for(j=0; j<NUM_DIM; j++)
				centroids[ pointBelongsToCluster[i] ][j] += points[i][j];

		// Calculate mean by dividing through number of points in each cluster
		for(k=0; k<NUM_CLUSTERS; k++) {
			for(j=0; j<NUM_DIM; j++)
				centroids[k][j] /= clusterContainsNumPoints[k];
		}

		/*
			2) Perform clustering using update strategy U-2: MacQueen.
		*/

		for(iter=0; iter<MAX_ITER; iter++) {

			membershipChange = false;

			for(i=0; i<NUM_POINTS; i++){
				// a) Reassign one point to another cluster (if there is a cluster center which is closer than the center of the current cluster)
				smallestSqDist = Float.MAX_VALUE;
				for(k=0; k<NUM_CLUSTERS; k++) {
					currentSqDist = 0.0f;
					for(j=0; j<NUM_DIM; j++) {
						diff = points[i][j] - centroids[k][j];
						currentSqDist += diff * diff;
					}

					if(currentSqDist < smallestSqDist) {
						smallestSqDist = currentSqDist;
						nearestCentroid = k;
					}
				}

				// b) Remember previous cluster membership of the point that is updated now, and update number of points in relevant clusters
				previouslyAssignedCluster = pointBelongsToCluster[i];
				newlyAssignedCluster = nearestCentroid;

				// check if point has moved from one to another cluster
				if(newlyAssignedCluster != previouslyAssignedCluster) {

					membershipChange = true;

					pointBelongsToCluster[i] = nearestCentroid;

					N_previouslyAssignedCluster = clusterContainsNumPoints[previouslyAssignedCluster];
					N_newlyAssignedCluster = clusterContainsNumPoints[newlyAssignedCluster];
					clusterContainsNumPoints[previouslyAssignedCluster]--;
					clusterContainsNumPoints[newlyAssignedCluster]++;

					// c) Determine new cluster center of every cluster in which a change has happened

					// Calculate the sum over all points in each cluster (but only of clusters with changes)
					for(j=0; j<NUM_DIM; j++) {
						centroids[newlyAssignedCluster][j] =
							( centroids[newlyAssignedCluster][j] * N_newlyAssignedCluster + points[i][j] )
								/ clusterContainsNumPoints[newlyAssignedCluster];
						centroids[previouslyAssignedCluster][j] =
							( centroids[previouslyAssignedCluster][j] * N_previouslyAssignedCluster - points[i][j] )
								/ clusterContainsNumPoints[previouslyAssignedCluster];
					}
				}
			}

			if(!membershipChange) {
				break;
			}
		}
		return iter;
	}


	/**
	 * Main
	 */
	public static void main(final String[] args) throws IOException {

		/*
			Variable declarations
		*/
		int iter, number_of_nonempty_buckets, d;
		long startTimeReading, endTimeReading, startTimeHashing, endTimeHashing, startTimeClustering, endTimeClustering, startTimeWriting, endTimeWriting;

		// declare & allocate & initialize (zeros by Java specs) arrays
		final float[][] points = new float[NUM_POINTS][NUM_DIM];
		final float[][] centroids = new float[NUM_CLUSTERS][NUM_DIM];
		final int[] pointBelongsToCluster = new int[NUM_POINTS];
		final int[] clusterContainsNumPoints = new int[NUM_CLUSTERS];

		float [][] bucket_points;
		int[] bucket_pointBelongsToCluster;

		/*
			I) Read data from provided CSV file into data structures
		*/

		startTimeReading = System.nanoTime();
		read_points_from_file( points );
		endTimeReading = System.nanoTime();

		/*
			II) Calculate a clustering with k-means
		*/

		startTimeHashing = System.nanoTime();

		// 0) Initialize hashes/buckets
		init_hashfunctions();
		number_of_nonempty_buckets = points2buckets( points );

		// prepare bucket_points for k-means, i.e. store copies of bucket
		// representatives in bucket_points
		bucket_points = new float[number_of_nonempty_buckets][NUM_DIM];
		assert number_of_nonempty_buckets == bucket2RepPointIdxMap.size();
		for( final Map.Entry<Integer,Integer> entry : bucket2bucketIdxMap.entrySet() ) {
			final Integer bucket = entry.getKey();
			final Integer bucketIdx = entry.getValue();
			assert bucketIdx < number_of_nonempty_buckets;
			// retrieve coordinates of representative
			final Integer pointIdx = bucket2RepPointIdxMap.get(bucket);
			for( d = 0; d < NUM_DIM; d++ )
				bucket_points[bucketIdx][d] = points[ pointIdx ][d];
		}

		endTimeHashing = System.nanoTime();

		startTimeClustering = System.nanoTime();

		// Initialize data structures that store information about our buckets
		bucket_pointBelongsToCluster = new int[number_of_nonempty_buckets];

		// 1) Perform clustering
		iter = kmeans( bucket_points, centroids, clusterContainsNumPoints, bucket_pointBelongsToCluster );

		// map bucket clustering results back to original points
		for( int p = 0; p < NUM_POINTS; p++ ) {
			final int bucketIdx = bucketIdxOfPoint[p];
			//System.out.println(String.format("point # %d -> bucket idx %d", p, bucketIdx));
			pointBelongsToCluster[p] = bucket_pointBelongsToCluster[ bucketIdx ];
		}
		endTimeClustering = System.nanoTime();

		/*
			III) Write data labeled with classes to CSV file
		*/

		startTimeWriting = System.nanoTime();
		write_results_to_file( pointBelongsToCluster );
		endTimeWriting = System.nanoTime();

		/*
			IV) Report timings and results.
		*/

		// print bucket analysis
		System.out.println( "Bucket analysis" );
		System.out.println( "---------------" );
		System.out.println(String.format("Number of non-empty buckets: %d [%.2g%%]", number_of_nonempty_buckets, 100f*number_of_nonempty_buckets/NUM_BUCKETS));
		print_clustering_analysis(points, bucketIdxOfPoint, "bucket");
		System.out.println();
		System.out.println("============================================================");

		System.out.println();
		System.out.println("Final clustering result");
		System.out.println("-----------------------");

		if(iter < MAX_ITER)
			System.out.println("Clustering converged in " + iter + " iterations.");
		else
			System.out.println("Clustering did not converge in a preset maximum of " + iter + " iterations.");

		print_clustering_analysis( points, pointBelongsToCluster, "cluster" );

		// Plot times
		System.out.println();
		System.out.println("============================================================");
		System.out.println();
		System.out.println("Time measurement");
		System.out.println("----------------");

		System.out.print("Elapsed time for loading data: ");
		System.out.println( ( endTimeReading-startTimeReading ) / 1e9f + "s");

		System.out.print("Elapsed time for hashing: ");
		System.out.println( ( endTimeHashing-startTimeHashing ) / 1e9f + "s");

		System.out.print("Elapsed time for clustering: ");
		System.out.println( ( endTimeClustering-startTimeClustering ) / 1e9f + "s");

		System.out.print("Elapsed time for writing data: ");
		System.out.println( ( endTimeWriting-startTimeWriting ) / 1e9f + "s");

		System.out.println();
		System.out.print("Elapsed time for hashing + clustering: ");
		System.out.println( ( endTimeClustering-startTimeHashing ) / 1e9f + "s");
	}
}
