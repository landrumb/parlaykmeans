The point of this experiment is to bench our different initialization functions. This benching has been done informally but this experiment will add graphs and csv files. 

Since Kmeans++ blows up with high values, 
we will run on 
k= 10, 100, 200, 400, 500, 1000, 2000, 3125, 4000, and 5000

to capture this blowup. Note that we will only run Kmeans++ on 10,100,200,400,500, and 1000.

In this way, this is an improvement on experiments 11 and 15.

Run on fern, August 8, 2023.

Run on this commit:
748f041e617bd69cb536777302e1c1a9f826243a