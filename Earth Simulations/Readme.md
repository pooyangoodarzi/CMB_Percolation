This is the simulation based on the paper https://arxiv.org/abs/1809.00256 done by Jingfang Fan, Jun Meng, Abbas Ali Saberi. The study is about exploring the elevations of the earth as a network, using percolation. The writers found out that the network has a remarkable feature which is explained by the theory of critical phenomena and phase transitions. they show that the network has an explosive percolation which I will explain a bit more later on. the idea is simple to assume the earth has no water on it then if I start to put water on it little by little and each time I calculate the biggest cluster of joint waters I can see that there exists a specific height where there exists the biggest possible cluster of waters. we can do this process in reverse as well if I fill the earth with water to the point that there exists no land and then little by little I take out the water and calculate the biggest cluster of land then again there exist a height where we can see the biggest land cluster. in this critical height if I take a little more water then there would be a lot of lands and if I put on a little more water I will lose many lands. the authors show that this critical height is the see level right now on earth. Me and my college (Pooyan) wrote the code together to drive the paper's results.
we used the data mentioned in the paper
the algorithm:
-1, sort all the elevation on earth in an ascending manner (for water clusters, if you want the land clusters sort the data in a descending manner)
-2, here each height is a node in our network
-3, we turn on each node in an increasing order
-4, we use the Hoshen-Kopelman algorithm with periodic boundary conditions on the zonal sides of the lattice to simulate the earth's spherical shape.
-5, we add the effect of the spherical form of the earth by calculating the cosine of the latitude of the used nodes in each step and we find the biggest cluster.
Since the network has explosive percolation, there is no finite-size effect for the different sizes of the lattice; this means that; I should have the same results even for big coursed grained lattices.
Here are the code results in comparison to the paper's results.

- Sort all the elevation on earth in an ascending manner (for water clusters, if you want the land clusters sort the data in a descending manner).
- Item 2
- Item 3

