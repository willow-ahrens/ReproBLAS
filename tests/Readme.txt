This directory contains a set of tests for the Reproducible Sum Algorithms:

Test 4 collects performance data for MPI-version
./run-test4
	-n size of vector
	-t mask of tests to be conducted:
		xxx1 Fast Reproducible Sum
		xx1x Blocked Reproducible Sum
		x1xx Normal dasum
		1xxx Accurate (Slow) Repdoucible Sum

Test 7 collects performance data running on 1 core for different vector size
for Normal dasum, Fast reproducible Sum, and Accurate (Slow) Reproducible Sum
./run-test7
	-n xx:xx:xx vector size

