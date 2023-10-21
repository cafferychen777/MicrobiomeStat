MicrobiomeStat: 支持R中的纵向微生物组分析
================

<p align="center" style="margin:0; padding:0;">
<img src="man/figures/logo.png" alt="MicrobiomeStat Logo" width="400" style="margin:0; padding:0;"/>
</p>
<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/MicrobiomeStat)](https://cran.r-project.org/package=MicrobiomeStat)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/MicrobiomeStat)](https://cran.r-project.org/package=MicrobiomeStat)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

[EN](README.md) \| [CN](README-CN.md) \| [ES](README-ES.md)

`MicrobiomeStat`包是一个专门用于探索**纵向微生物组数据**的R工具。它也适用于多组学数据和横断面研究，重视社区内的集体努力。此工具旨在支持研究人员在时间上进行广泛的生物学查询，对社区现有资源表示感谢，并以合作精神推动微生物组研究的进一步发展。

# 新闻

📢 **更新** (10月20日):
Shiny界面现已正式上线，目前支持小型至中型数据集的分析。您可以通过[此链接](https://microbiomestat.shinyapps.io/MicrobiomeStat-Shiny/)访问界面。

如果服务器限制影响了您的分析，或者您更倾向于使用较小的模块，我们建议您直接使用我们的包。对于需要更灵活部署的情况，可以考虑从[此处](https://microbiomestat.shinyapps.io/MicrobiomeStat-Shiny/)克隆我们的Shiny仓库，并将其部署在您的本地服务器或计算机上。

我们感谢您的理解和持续的参与。

# 引文

## 对MicrobiomeStat的一般引用

如果您使用的功能超出了`linda`和`linda.plot`函数，请按照以下方式引用，直到预印本版本发布：

    @Manual{,
      title = {MicrobiomeStat: Comprehensive Statistical and Visualization Methods for Microbiome and Multi-Omics Data},
      author = {Xianyang Zhang and Jun Chen and Caffery(Chen) Yang},
      year = {2023},
      note = {R package version 1.1.1},
      url = {https://www.microbiomestat.wiki},
    }

## 对特定`MicrobiomeStat`函数的引用

如果您正在使用`linda`，`linda.plot`，`generate_taxa_association_test_long`，`generate_taxa_test_pair`，`generate_taxa_test_single`或`generate_taxa_trend_test_long`函数，请引用以下论文：

    @article{zhou2022linda,
      title={LinDA: linear models for differential abundance analysis of microbiome compositional data},
      author={Zhou, Huijuan and He, Kejun and Chen, Jun and Zhang, Xianyang},
      journal={Genome biology},
      volume={23},
      number={1},
      pages={1--23},
      year={2022},
      publisher={BioMed Central}
    }

我们将在预印本发布后更新引用指南。

## 关于CRAN版本的重要说明

`MicrobiomeStat`包正在不断开发中。因此，最新的功能尚未纳入CRAN存储库中可用的版本。当前的CRAN版本仅支持`linda`和`linda.plot`函数。对于需要更广泛功能的用户，特别是与纵向数据分析相关的功能，建议直接从GitHub安装开发版本。这个过程需要先安装`devtools`包。

``` r
install.packages("devtools")
```

安装了`devtools`之后，您可以使用以下命令从GitHub安装`MicrobiomeStat`：

``` r
devtools::install_github("cafferychen777/MicrobiomeStat")
```

# 目录

1.  [引文](#citations)
    - [对MicrobiomeStat的一般引用](#general-citation)
    - [对特定`MicrobiomeStat`函数的引用](#specialized-citation)
    - [关于CRAN版本的重要说明](#cran-version-note)
2.  [在线教程](#online-tutorials)
    - [熟悉流程，享受无缝体验](#acquaint-yourself)
    - [探索MicrobiomeStat](#explore-microbiomestat)
      - [功能探索](#feature-exploration)
      - [致谢](#acknowledgements)
3.  [使用MicrobiomeStat的好处](#why-choose-microbiomestat)
    - [用户支持](#user-support)
    - [持续开发](#ongoing-development)
    - [合作开发](#collaborative-development)
    - [结论和特性](#conclusion)
    - [演示报告](#demo-reports)
4.  [帮助和联系信息](#support-contact)
5.  [参与我们的Discord社区](#discord-community)
6.  [分享和连接](#share-and-connect)

# 在线教程

`MicrobiomeStat`提供了一套全面的微生物组数据分析工具，包括从数据输入到可视化的各种功能。

为了让用户熟悉`MicrobiomeStat`，我们在GitBook上提供了详细的在线教程。教程涵盖了以下内容：

- 安装和配置说明
  - 这些指南有助于确保您的设置正确配置和优化。
- 基于现实场景的分析演示
  - 这些演示提供实用的见解和技能。
- 练习的代码示例
  - 这些示例让用户熟悉`MicrobiomeStat`的编码实践。
- 解读结果和创建可视化的指南
  - 这些指南帮助用户理解和有效地展示他们的数据。
- 常见问题的答案
  - 本节提供了对常见问题的快速解决方案。

## 熟悉流程，享受无缝体验

为了在`MicrobiomeStat`中获得无缝体验，最大限度地利用这些丰富的资源：

[**📘 探索MicrobiomeStat教程**](https://www.microbiomestat.wiki)

## 探索MicrobiomeStat

微生物组研究领域复杂且不断发展。选择的分析工具在研究过程中可以起到关键作用。在这种情况下，`MicrobiomeStat`旨在成为一种支持伙伴。

### 功能探索

为了更好地理解`MicrobiomeStat`的能力，我们在我们的网站上概述了其特性和功能，并提供了一些其他工具的上下文引用，以便有更多的信息：

- [探索纵向分析特性](https://www.microbiomestat.wiki/introduction/unveiling-microbiomestat-a-glimpse-into-diverse-microbiome-analysis-solutions/exploring-longitudinal-analysis-features-microbiomestat-among-other-packages)

- [探索集成分析特性](https://www.microbiomestat.wiki/introduction/unveiling-microbiomestat-a-glimpse-into-diverse-microbiome-analysis-solutions/exploring-feature-offerings-microbiomestat-among-integrated-analysis-packages)

### 致谢

我们在`MicrobiomeStat`上站在巨人的肩膀上，我们对依赖于我们包的勤奋和杰出的开发人员表示衷心的感谢。他们的卓越努力不仅使我们的工作成为可能，而且还显著提高了科学社区可用的计算工具的标准：

- 核心依赖：
  - R (\>= 3.5.0), rlang, tibble
- 导入的包：
  - ggplot2, matrixStats, lmerTest, foreach, modeest, vegan, dplyr,
    pheatmap, tidyr, ggh4x, ape, GUniFrac, scales, stringr, rmarkdown,
    knitr, pander, tinytex
- 建议的包：
  - ggrepel, parallel, ggprism, aplot, philentropy, forcats, yaml,
    biomformat, Biostrings

此外，我们对在微生物组研究社区中创建和维护以下卓越工具的开创者表示最深的赞赏和尊重。他们的开创性工作已经在微生物组数据分析的复杂景观中铺设了路径，我们真心荣幸能够与他们并肩前行：

- `microbiomeutilities`, `phyloseq`, `microbiomemarker`,
  `MicrobiomeAnalyst`, `microbiomeeco`, `EasyAmplicon`, `STAMP`,
  `qiime2`, 和 `MicrobiotaProcess`

他们的贡献激励我们继续改进和扩展`MicrobiomeStat`的功能，我们真诚地希望我们的谦卑添加能够成为已经为研究者提供的令人难以置信的工具阵列的有用补充。

### 用户支持

`MicrobiomeStat`设计时考虑到用户。我们提供全面的[文档和教程](https://www.microbiomestat.wiki/)，以协助新手和有经验的研究者。在发布问题或问题之前，我们鼓励用户[检查以前的问题和问题](https://github.com/cafferychen777/MicrobiomeStat/issues?q=is%3Aissue+is%3Aclosed)，看看该主题是否已经得到了解答。

如果您对特定函数的文档有特定的评论或问题，而您发现RStudio的搜索框导致404错误，您可以直接在<https://cafferychen777.github.io/MicrobiomeStat/reference/index.html>访问函数的文档。

如果您的问题或疑问尚未得到解答，随时在GitHub上[开启新问题](https://github.com/cafferychen777/MicrobiomeStat/issues)。我们非常鼓励全球的研究者、开发者和用户之间的交流与合作。为了促进全球间的沟通和互助，我们希望您能使用英文提问。这样不仅可以帮助更多的人理解和参与到讨论中，还能确保我们更迅速、准确地为您提供帮助。谢谢您的理解和支持！我们在这里，随时为您解决可能遇到的任何挑战。

### 持续开发

确保`MicrobiomeStat`保持其类别中的领先地位需要持续的开发。我们致力于定期更新和解决用户反馈。作为我们持续改进努力的一部分，我们为`MicrobiomeStat`开发了一个Shiny界面，它提供了一种用户友好、交互式的方式来进行微生物组数据分析。Shiny界面与主包一起得到了积极的维护和改进。你可以在其专用的[GitHub仓库](https://github.com/cafferychen777/MicrobiomeStat-Shiny)中访问Shiny应用程序文件和本地设置指南。

### 合作开发

`MicrobiomeStat`是一个开源工具，我们非常重视来自社区的贡献。如果你有建议，改进或对未来开发方向和功能添加的反馈，[拉取请求](https://github.com/cafferychen777/MicrobiomeStat/pulls)是受欢迎的，你也可以在我们的GitHub仓库的[讨论区](https://github.com/cafferychen777/MicrobiomeStat/discussions)分享你的想法。与其他社区成员互动，帮助我们使`MicrobiomeStat`成为微生物组研究的更有用的工具。

### 结论

`MicrobiomeStat`旨在为微生物组数据分析提供一种可靠和高效的资源。我们诚邀所有重视开源合作的人加入我们的社区，共同为其持续发展做出贡献。

| 功能                                                                                                                                                                        | 描述                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [数据导入和转换](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object) | 支持从[QIIME2](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object/importing-data-from-qiime2-into-microbiomestat)，[Mothur](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object/fetching-data-from-mothur-into-microbiomestat)，[DADA2](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object/integrating-data-from-dada2-into-microbiomestat)，[Phyloseq](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object/navigating-data-from-phyloseq-into-microbiomestat)等平台的多种输入格式。 |
| [横断面研究分析](https://www.microbiomestat.wiki/cross-sectional-study-design/unraveling-cross-sectional-studies-with-microbiomestat)                                       | 为横断面研究提供全面的分析。                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| [配对样本分析](https://www.microbiomestat.wiki/paired-samples-analysis/unveiling-paired-samples-analysis-a-comprehensive-guide)                                             | 提供配对样本分析的工具。                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| [纵向研究分析](https://www.microbiomestat.wiki/longitudinal-study-design/grasping-longitudinal-studies-introduction-and-dataset-overview)                                   | 有助于探索微生物组的时间动态。                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| 报告生成功能                                                                                                                                                                | 包括为[横断面](https://www.microbiomestat.wiki/cross-sectional-study-design/cross-sectional-reporting-microbial-analysis-reports-with-microbiomestat)，[配对](https://www.microbiomestat.wiki/paired-samples-analysis/automated-reporting-for-paired-studies-microbiomestats-integrated-analysis-reports)，[纵向](https://www.microbiomestat.wiki/longitudinal-study-design/longitudinal-reporting-microbiome-analysis-automation-with-microbiomestat)研究设计提供单独的报告功能。正在开发Shiny接口，以实现一键生成报告。                                                                                                                                                                                                                                                                                                                                                                                      |
| 可视化能力                                                                                                                                                                  | 支持广泛的可视化样式。                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| 持续开发                                                                                                                                                                    | 致力于持续完善现有功能并添加新功能。                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |

这种方法保证了用户可以轻松导航到`MicrobiomeStat`文档的特定部分，获取各种分析类型的详细信息和指南。这种结构和可访问性帮助用户有效地利用`MicrobiomeStat`满足他们的微生物组数据分析需求。

### 演示报告

对于那些有兴趣看到`MicrobiomeStat`实际操作的人，我们已经准备了针对不同研究设计的演示报告：

- [横断面研究设计：使用MicrobiomeStat报告微生物分析](https://www.microbiomestat.wiki/cross-sectional-study-design/cross-sectional-reporting-microbial-analysis-reports-with-microbiomestat)
- [配对样本分析：使用MicrobiomeStat报告微生物分析](https://www.microbiomestat.wiki/paired-samples-analysis/automated-reporting-for-paired-studies-microbiomestats-integrated-analysis-reports)
- [纵向研究设计：使用MicrobiomeStat进行微生物分析自动化](https://www.microbiomestat.wiki/longitudinal-study-design/longitudinal-reporting-microbiome-analysis-automation-with-microbiomestat)

我们鼓励你探索这些示例，并发现我们工具的强大功能。

## 帮助和联系信息

如需帮助或查询，欢迎联系：

| 姓名                                                                         | 邮箱                        |
|------------------------------------------------------------------------------|-----------------------------|
| [Dr. Jun Chen](https://scholar.google.com/citations?user=gonDvdwAAAAJ&hl=en) | <Chen.Jun2@mayo.edu>        |
| [Chen YANG](https://cafferyyang.owlstown.net/)                               | <cafferychen7850@gmail.com> |

## 加入我们的Discord社区

加入我们的Discord社区，了解`MicrobiomeStat`的最新更新、发展和增强。参与热烈的讨论，提问，分享见解，从同行和专家那里获得支持：

[加入MicrobiomeStat Discord服务器！](https://discord.gg/BfNvTJAt)

在我们的Discord服务器中，一个自动化的机器人会将每个包和教程更新的信息告知你，确保你不会错过新功能、改进和学习材料。我们活跃的社区依赖于合作、反馈和持续学习，对于初学者和经验丰富的研究人员在微生物组数据分析世界中导航都是无价的空间。保持联系，保持了解，让我们一起推进微生物组数据分析领域的发展！

## 分享和连接

通过各种平台传播关于`MicrobiomeStat`的信息并保持连接！

[![Twitter](https://img.shields.io/twitter/url?url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&style=social)](https://twitter.com/intent/tweet?url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&text=Check%20out%20this%20awesome%20package%20for%20comprehensive%20and%20longitudinal%20microbiome%20analysis%20in%20R!)

[![Facebook](https://img.shields.io/badge/Share_on-Facebook-1877F2?logo=facebook&style=social)](https://www.facebook.com/sharer/sharer.php?u=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&quote=Check%20out%20this%20awesome%20package%20for%20comprehensive%20and%20longitudinal%20microbiome%20analysis%20in%20R!)

[![LinkedIn](https://img.shields.io/badge/Share_on-LinkedIn-0077B5?logo=linkedin&style=social)](https://www.linkedin.com/shareArticle?mini=true&url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&title=Check%20out%20this%20awesome%20package%20for%20comprehensive%20and%20longitudinal%20microbiome%20analysis%20in%20R!)

[![Reddit](https://img.shields.io/badge/Share_on-Reddit-FF4500?logo=reddit&style=social)](https://www.reddit.com/submit?url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&title=Check%20out%20this%20awesome%20package%20for%20comprehensive%20and%20longitudinal%20microbiome%20analysis%20in%20R!)

[![WhatsApp](https://img.shields.io/badge/Share_on-WhatsApp-25D366?logo=whatsapp&style=social)](https://wa.me/?text=Check%20out%20this%20awesome%20package%20for%20comprehensive%20and%20longitudinal%20microbiome%20analysis%20in%20R!%20https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat)

[![Slack](https://img.shields.io/badge/Share_on-Slack-4A154B?logo=slack&style=social)](https://slack.com/intl/en-cn/)

[![Email](https://img.shields.io/badge/Share_on-Gmail-D14836?logo=gmail&style=social)](mailto:?subject=Check%20out%20this%20awesome%20package%20for%20comprehensive%20and%20longitudinal%20microbiome%20analysis%20in%20R!&body=I%20found%20this%20amazing%20R%20package%20for%20microbiome%20analysis%20called%20%60MicrobiomeStat%60.%20Check%20it%20out%20here:%20https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat)
