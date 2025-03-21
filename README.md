# Effects of formalin fixation on polarimetric properties of brain tissue: fresh or fixed?

by
Éléa Gros, Omar Rodrıguez-Nunez, Leonard Felger, Stefano Moriconi, Richard McKinley, Angelo Pierangelo, Tatiana Novikova, Erik Vassella, Philippe Schucht, Ekkehard Hewer, Theoni Maragkou

This paper has been accepted for publication in *Neurophotonics*.

This study aimed to analyse the effects of formalin fixation on the polarimetric parameters of brain tissue to evaluate the probability of success of the future transfer of ML algorithms from fixed to fresh brain tissue and vice-versa.


## Abstract

**Significance**: Imaging Mueller polarimetry (IMP) appears as a promising technique for real-time delineation of healthy and neoplastic tissue during neurosurgery. The training of machine learning algorithms used for the image post-processing requires large data sets that can be derived from the measurements of formalin-fixed brain sections. However, the success of the transfer of such algorithms from fixed to fresh brain tissue and vice-versa depends on the degree of alterations of polarimetric properties induced by formalin fixation (FF).

**Aim**: The comprehensive studies were performed on the FF induced changes in fresh brain tissue polarimetric properties.

**Approach**: Polarimetric properties of pig brain were assessed in thirty coronal thick sections before and after FF using a wide-field IMP system. The width of the uncertainty region between grey and white matter was also estimated. 

**Results**: The depolarization increased by 2% in grey matter and decreased by 1% in white matter following FF, whereas linear retardance decreased by 30% in both tissue types after FF. The visual contrast between grey and white matter and fiber tracking remained preserved after FF. Tissue shrinkage induced by FF did not have significant effect on the uncertainty region width.

**Conclusions**: Similar polarimetric properties were observed in both fresh and fixed brain tissues, indicating a high potential for transfer learning.


## Software implementation

This GitHub folder documents all the source code used to generate the results and figures in the paper.

All source code used to generate the results and figures in the paper are in the `base` folder, in three different [Jupyter notebooks](http://jupyter.org/): `selection_of_ROIs.ipynb`,  `show_histograms.ipynb`, `parameter_comparaison.ipynb` and `evaluation_border_zone.ipynb`. The data used in this study should be copied from [this link](to be changed), and placed in the `data` folder. Results generated by the code are saved in `results`.


## Getting the code

You can download a copy of all the files in this repository by cloning the
[git](https://git-scm.com/) repository:

    git clone https://github.com/RomGr/FixationEffectPaper.git


## Selection of ROIs

The first juypter notebook can be run to select automatically 25 ROIs within the white matter, and 25 ROIs for each series of measurement located in the folder `data/fixation_over_time`. CAUTION: running the notebook erase the previously generated results.


## Show histograms

The second juypter notebook allows to visualize the distributions of several ROIs, and to determine wether or not these distributions are normal. The output for this notebook can be found in `results/histograms`.


## Parameter comparison

The third juypter notebook, which should be run after the selection of the ROIs, was used to perform the statistical analyses for the evolution of the means and the fold changes. The results generated (mainly excel tables), can be found in the folder `results/comparaison`. The plots for the manuscript were created using Graphpad Prism.


## Evalutation border zone

The last juypter notebook, independent for the first two ones, allows to reproduce the results for the evalutation of the uncertainty region in between white and grey matter, and its evolution following formalin fixation. The results generated can be found in the folders `results/fixed` and `results/fresh`. The plots for the manuscript were created using python and are available directly in the `results` folders.


## Data

Three subfolders can be found in the `data` folder:
1. `fixation_over_time`: contains the measurements for the section performed at different time points
2. `fresh`: contains the measurements for the fresh section
3. `fixed`: contains the measurements for the fixed section performed 24 hours after the fresh ones


## License

All source code is made available under a MIT license. See `LICENSE` for the full license text.
