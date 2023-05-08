![ExP-MRM标志](http://www.exposomemrm.com/static/img/website_icon.dfecee9.png "CIL-MRM logo")

<!-- TOC -->

- [ExP-MRM introduction](#ExPMRM introduction)
    - [predictionRT.R](#predictionRT.R)
    - [transitiongroupdata.R](#transitiongroupdata.R)
    - [derivatization.py](#derivatization.py)
    - [peakanalysis.R](#peakanalysis.R)
    - [CIL-PMRM Exposome Database](#CIL-PMRM Exposome Database)
- [install](#install)
    - [ExP-MRM OP](#ExP-MRM OP)
    - [ExP-MRM ENV](#ExP-MRM ENV)
- [快速入门](#快速入门)
    - [predictionRT](#predictionRT)
    - [transitiongroupdata](#transitiongroupdata)
    - [peakanalysis](#peakanalysis)
- [联系我们](#联系我们)
- [官方](#官网)
- [许可证](#许可证)

<!-- /TOC -->

## ExPMRM introduction

The ExPMRM project is built using Python and R, leveraging a variety of libraries. For the backend in silico derivatization of targeted compounds, the RDKit Python library is employed. It is important to mention that the derivatization module has been developed as a standalone application within this project. The backend computations for retention time (RT) prediction and mass spectrometry (MS) data analysis rely on several R packages, including the rcdk package, the randomForest package, and the xmcs package. 


### predictionRT.R

We adopted a RT prediction model based on random forest algorithm, by using available standards to obtain the RT windows of 110 k compounds in the database. The predicted RT window of DnsCl-OH derivatization products was 1.26 min with R = 0.85 for the total running time of 20 min. Similarly, the RT window of 174 MPEA-COOH derivatization products was 0.91 min with R = 0.93.


### transitiongroupdata.R
We generated a dynamic MRM optimization algorithm to ensure efficient detection. Compounds with the same m/z and the difference value of their m/z is 2 Da (OH: m/z of 13C2-DnsCl — m/z of DnsCl), or 3 Da (COOH: m/z of d3-MPEA — m/z of MPEA) were randomly assigned to different CIL-pseudo-MRM methods as the format of csv. MRM transition capacity was another issue should be taken into consideration. For AB QTRAP 6500, the dwell time at 3 ms was enough, so each segment could contain almost 250 MRM transitions.

### derivatization.py
The chemical structures of MPEA-COOH and DnsCl-OH derivatized products were generated using RDKit. Subsequently, we searched in their SMILES structures to observe whether there were any -COOH or -OH groups using SMARTS pattern matching. In that case, the -OH group was replaced by the SMILES structure of DnsCl (CN(C)C1=CC=CC2=C([S](=O)=O)C=CC=C21); while the -COOH group was replaced by the SMILES structure of MPEA (C[N]CCC1=CC=CC=C1).

### peakanalysis.R
MS data were uploaded after LC-MS acquisition. Compounds with the ΔRT of light- and heavy-derivatized MRM transitions < 0.1 min and intensity ratios within 0.5 – 2.0, were extracted and considered as potential OH or COOH compounds. Then, t-test and ANOVA were used for the data statistical difference analysis of peak intensity of the potential environmental biomarkers in biological samples. Compounds with p value < 0.05 and fold change > 1.5 were defined as significantly changed compounds.

### CIL-PMRM Exposome Database
The ExPMRM database contains chemicals with great environmental concern, large production, high human exposure or toxicity. In the database, 2612 parent compounds containing OH or COOH were acquired from HExpMetDB and their > 110 k biotransformation metabolites were predicted by BioTransformer 1.04. In this study, we have developed a ExPMRM database of > 110 k compounds with OH or COOH group by CIL strategy. The database contained several items of environmental pollutants, such as name, InChIKey, RT, SMILES structure, exact precursor m/z, product ion, and CE values, some of which were obtained by experimental results and related R and Python packages.

## install

### ExP-MRM OP

|  OP       | status  |
|  :-------------- | :--- |
| windows     | ✔️   |
| ubuntu-x86  | ✔️   |


### ExP-MRM ENV
|  language        | version  |
|  :-------------- | :--- |
| python    | 3.6+  |
| R  |  4.2.0  |

1. Please download and install the WHL package from the official installation.

R language package
|  package        | version  |
|  :-------------- | :--- |
| rcdk    | 4.1.1  |
| rcdklibs  |  4.0.0 |
| dplyr    | 4.2.0 |
| rio  |  4.0.0  |
| hash    | 4.2.0  |
| xcms  |  4.2.0  |
| magrittr  |  4.1.3  |
| MSnbase    | 4.2.0  |
| tidyr  |  4.1.2  |
| ggplot2  |  4.2.0  |
| tidyverse    | 4.2.0   |
| ggpubr  |  4.2.0  |
| ggrepel  |  4.2.0  |
| ggfortify    | 4.1.3  |

python language package
|  package        | version  |
|  :-------------- | :--- |
| kora    | 0.9.20  |
| pandas  |  1.4.3+ |


2. [建议使用anaconda环境按照相关rdk依赖包](https://github.com/rdkit/rdkit/blob/master/Docs/Book/Install.md)

3. [ubuntu建议使用apt-get安装R语言依赖包](http://ftp.sjtu.edu.cn/ubuntu/pool/universe/r)


4. 执行以下命令，验证安装结果。

    ```python
    import kora as kora
    import pandas as pd
    ```

    ```R
    packsneed <- c('xcms','magrittr','MSnbase','dplyr','tidyr','ggplot2','tidyverse','ggpubr',"ggrepel","rio", 'ggfortify')
    packsneed <- c('rcdk','rcdklibs','dplyr',"rio", 'hash')
    ```


## install

参考[install](http://www.exposomemrm.com/about)实现。
### predictionRT
    ```R
    Rscript ./predictionRT/predictionRT.R output.csv OH ./output
    ```
### transitiongroupdata
    ```R
    Rscript transitiongroup.R input.csv 50 ./output
    ```
### peakanalysis
    ```R
    Rscript peakanalysis.R  infilename1 infilename2   0.01 30 0.1 parameter 300 ./output
    ```


## 官网
参考[官网](http://www.exposomemrm.com)

## 联系我们
有关安装指南、教程和API的更多详细信息，请联系我们(http://www.exposomemrm.com/contact)。

## 许可证
[Apache License 2.0](https://gitee.com/mindspore/mindspore/blob/master/LICENSE)
