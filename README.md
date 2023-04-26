![CIL-MRM标志](http://www.exposomemrm.com/static/img/website_icon.dfecee9.png "CIL-MRM logo")

<!-- TOC -->

- [CIL-MRM介绍](#CIL-MRM介绍)
    - [Pseudo-multiple reaction monitoring](#Pseudo-MRM)
    - [Derivatization-LC-MS strategy](#Derivatization-LC-MS)

- [安装](#安装)
- [快速入门](#快速入门)
- [功能介绍](#功能介绍)
   - [peakanalysis](#peakanalysis)
   - [peakanalysis](#peakanalysis)
   - [peakanalysis](#peakanalysis)
- [联系我们](#联系我们)
- [官方](#官网)
- [许可证](#许可证)

<!-- /TOC -->

## CIL-MRM介绍

An automatic and integrative exposome platform
>
欲了解更多详情，请查看我们的[官网](http://www.exposomemrm.com)。

### Pseudo-MRM

It is usually developed by transforming anuntargeted compounds profiling method to a pseudo-targeted MRM method. In detail, the ion pairs of pseudo-MRM of these methods were acquired from the real samples through either untargeted tandem MS/MS by HRMS or directly from MS/MS compound database.Pseudo-multiple reaction monitoring (Pseudo-MRM) merges the advantages of both untargeted (high coverage) and targeted detection (high sensitivity) methods.

## Derivatization-LC-MS
Derivatization is also named chemical isotope labeling. In such a strategy, a pair of isotope-codedreagents are used to tag with reactive functional groups before generated products are subjected to LC-MS analysis.The strategy can further increase the detection sensitivity of analytes, improve their chromatographic peak shapes and provide one-to-one internal standards to reduce false positive rates.


## 安装

CIL-MRM提供跨多个后端的构建选项：

|  操作系统        | 状态  |
|  :-------------- | :--- |
| windows 系列    | ✔️   |
| ubuntu-x86  | ✔️   |


CIL-MRM基础环境设置：
|  安装包        | 版本  |
|  :-------------- | :--- |
| python    | 3.6+  |
| R  |  4.2.0  |

1. 请从官方安装下载并安装whl包。

R语言包依赖
|  安装包        | 版本  |
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

python 语言依赖包
|  安装包        | 版本  |
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


## 快速入门

参考[快速入门](http://www.exposomemrm.com/about)实现。

## 功能介绍
### peakanalysis

## 官网
参考[官网](http://www.exposomemrm.com)

## 联系我们
有关安装指南、教程和API的更多详细信息，请联系我们(http://www.exposomemrm.com/contact)。



## 许可证

[Apache License 2.0](https://gitee.com/mindspore/mindspore/blob/master/LICENSE)
