###########################
*tempowork* Python Library
###########################

--------------------------------------------------
 Frequent subgraph mining in temporal networks
--------------------------------------------------

The *tempowork* Python library comprises a series of algorithms for mining temporal networks where time is represented as a continuous dimension. In the current version, given a user-defined frequency threshold ``supp`` in  ``{1,2,...,n}``, it can mine all the frequent subgraphs in a data set of ``n`` temporal networks. The frequent subgraphs are subgraphs that appear in at least ``supp`` networks of the data set. The mining can be executed under one of the four following definitions of isomorphism: 

	* Exact-time isomorphism
	* Inexact-time isomorphism
	* Sequence-preserved exact-time isomorphism
	* Sequence-preserved inexact-time isomorphism

--------------------------------------------------
 Mining evolving patterns in temporal networks
--------------------------------------------------
Coming soon! 

*************
Install
*************

Install the latest version of tempowork:

.. code-block:: python

	pip install tempowork

*************
Usage
*************

The main function can be executed as follows:

.. code-block:: python

	import tempowork as tw
	test = tw.frequent(file_name, support, isomorphism = 'e', disc_threshold = 0.05, binning=False, nbins = 10, directed = False)

**Arguments**

* *file_name*: a string referring to the data set of temporal networks (`.txt` file). 
* *support*: the support integer value in ``{1,2,...,n}``, representing the number of networks in the data set that should include a pattern to make the pattern frequent
* *isomorphism*: it shows when two networks are considered isomorphic. It can take one of the following options: 

  * *'e'*: Exact-time isomorphism
  * *'i'*: Inexact-time isomorphism
  * *'es'*: Sequence-preserved exact-time isomorphism
  * *'is'*: Sequence-preserved inexact-time isomorphism

* *disc_threshold*: a user-defined threshold used for 'i' and 'is' isomorphism definitions. Once all the durations are sorted, disc_threshold is the maximum value that the difference between two consecutive durations divided by the first duration is permitted, in order to consider two consecutive values inexactly identical.
* *binning*: used for 'i' and 'is' isomorphism definitions and takes one of the following values.

  * *False*: the inexact identical durations are determined by disc_threshold.
  * *width*: the algorithm uses equal-width periods to find inexact identical values
  * *frequency*:  the algorithm uses an equal-frequency strategy to find inexact identical values

* *nbins*: it is a parameter that works with 'i' and 'is' isomorphism definitions when binning is not False. It determines the number of bins for 'width' and 'frequency' options of binning.
* *directed*: it specifies whether the (edges of) networks in the data set are directed. 


**Important**

	Each network in the temporal network data set comprises a list of temporal edges, ordered based on their starting time. Each line in the data set represents either a sequential number for a network, ``t # net_id``, or a temporal edge, ``e v1_id v2_id v1_lbl e_lbl v2_lbl st dt``, where ``st`` and ``dt`` are starting point and duration of the edge, respectively. Note that all the identifiers, labels, time points, and durations should be integer values. The following few lines show an example of a data set of temporal networks composed of two networks, each consisting of four temporal edges.
   
		.. code-block:: python
		
			t # 0
			e 0 1 1 1 1 0 80
			e 1 2 1 3 1 20 70
			e 0 2 1 2 1 30 60
			e 1 3 1 4 1 60 40
			t # 1
			e 0 1 1 1 1 0 80
			e 1 2 1 3 1 20 70
			e 0 2 1 2 1 30 60
			e 1 3 1 4 1 60 40
		
*************
Examples
*************

Here are a few examples of using the `tempowork` library to mine frequent subgraphs by adopting different isomorphism definitions and parameters. The corresponding text files are provided in the ``test`` folder.

.. code-block:: python

	import tempowork as tw
	exact_example = tw.frequent('exact.txt', 2, isomorphism = 'e')
	inexact_example = tw.frequent('inexact.txt', 2, isomorphism = 'i', disc_threshold = 0.05)
	seq_exact_example = tw.frequent('seq_exact.txt', 2, isomorphism = 'es')
	seq_inexact_example = tw.frequent('seq_inexact.txt', 2, isomorphism = 'is', disc_threshold = 0.5)
	seq_inexact_example = tw.frequent('seq_inexact.txt', 2, isomorphism = 'is', binning = 'width', nbins = 10)
	seq_inexact_example = tw.frequent('seq_inexact.txt', 2, isomorphism = 'is', binning = 'frequency', nbins = 10)


Then, the results can be examined using:

.. code-block:: python

	number_of_frequent_patterns = exact_example.frequent_cntr
	frequent_patterns_detected = exact_example.frequent_patterns


*******************************************************
Request for feedback (It remains a work in progress!)
*******************************************************

The implementation of this algorithm requires multiple components, such as interval trees, constrained interval graphs, different definitions of isomorphism, and ..., to work seamlessly together. I tried to implement them accordingly. So, if you encounter any strange behavior, I would be happy to hear about your experience for further improvements. Please feel free to reach out via email (ali.jazayeri@drexel.edu).


*************
Related work
*************
Some of the other algorithms in the literature are surveyed in the following two papers:

If the data set composed of a set or sequence of static or temporal networks (**Note:** this paper comes with supplementary materials):

	Jazayeri A. and Yang C. C., *Frequent Subgraph Mining Algorithms in Static and Temporal Graph-Transaction Settings: A Survey*, in IEEE Transactions on Big Data, 2021
	https://doi.org/10.1109/TBDATA.2021.3072001

If the data set represents one single large static or temporal network:
	
	Jazayeri A. and Yang C. C., *Motif Discovery Algorithms in Static and Temporal Networks: A Survey*, Journal of Complex Networks, Volume 8, Issue 4, 2020, cnaa031
	https://doi.org/10.1093/comnet/cnaa031

If you could not access these papers, please contact me.

*************
Citation
*************

**Paper:** To Be Provided!
