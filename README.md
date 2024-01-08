MicrobiomeStat: Supporting Longitudinal Microbiome Analysis in R
================

<p align="center" style="margin:0; padding:0;">
<img src="man/figures/logo.png" alt="MicrobiomeStat Logo" width="400" style="margin:0; padding:0;"/>
</p>

🌟 **If you find `MicrobiomeStat` helpful, please consider giving us a
star on GitHub!** Your support greatly motivates us to improve and
maintain this project. 🌟

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/MicrobiomeStat)](https://cran.r-project.org/package=MicrobiomeStat)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/MicrobiomeStat)](https://cran.r-project.org/package=MicrobiomeStat)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

[EN](README.md) \| [CN](README-CN.md) \| [ES](README-ES.md)

The `MicrobiomeStat` package is a dedicated R tool for exploring
**longitudinal microbiome data**. It also accommodates multi-omics data
and cross-sectional studies, valuing the collective efforts within the
community. This tool aims to support researchers through their extensive
biological inquiries over time, with a spirit of gratitude towards the
community’s existing resources and a collaborative ethos for furthering
microbiome research.

# News

📢 **Update** (January 8th): Enhancement in Color Palette Functionality

We are excited to announce a significant update to our color palette
management in the majority of our functions. Users can now directly use
predefined palette names such as “lancet”, “nejm”, among others, as
values for the `palette` parameter. This update is part of our ongoing
effort to enhance user experience and provide more intuitive and
flexible options for data visualization.

Key highlights of this update: - Integration of the `mStat_get_palette`
function across various functions in our package. - Allows users to
easily specify color palettes by name, such as “lancet”, “nejm”, “npg”,
“aaas”, “jama”, “jco”, and “ucscgb”. - Ensures backward compatibility
and introduces a more user-friendly approach to selecting color schemes
for visualizations. - Aims to enhance the visual appeal and
interpretability of plots and heatmaps generated using our package.

This update is immediately available and we encourage users to explore
these new options in their analyses. For detailed usage, please refer to
the updated function documentation.

📢 **Update** (October 20th): The Shiny interface is now officially
available for use. It is currently configured to handle analysis for
small to medium-sized datasets. The interface can be accessed via [this
link](https://microbiomestat.shinyapps.io/MicrobiomeStat-Shiny/).

In the event of server limitations affecting your analysis, or for those
preferring to work with smaller modules, we recommend using our package
directly. For a more flexible deployment, consider cloning our Shiny
repository from
[here](https://github.com/cafferychen777/MicrobiomeStat-Shiny) and
deploying it on your local server or computer.

We appreciate your understanding and continued engagement.

# Citations

## General Citation for MicrobiomeStat

If you are using features beyond the `linda` and `linda.plot` functions,
please cite as follows, until a preprint version is published:

    @Manual{,
      title = {MicrobiomeStat: Comprehensive Statistical and Visualization Methods for Microbiome and Multi-Omics Data},
      author = {Xianyang Zhang and Jun Chen and Caffery(Chen) Yang},
      year = {2023},
      note = {R package version 1.1.1},
      url = {https://www.microbiomestat.wiki},
    }

## Citation for Specialized `MicrobiomeStat` Functions

If you are using the `linda`, `linda.plot`,
`generate_taxa_association_test_long`, `generate_taxa_test_pair`,
`generate_taxa_test_single`, or `generate_taxa_trend_test_long`
functions, please cite the following paper:

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

We will update the citation guidelines as soon as the preprint is
published.

## Important Note on CRAN Version

The `MicrobiomeStat` package is under continuous development. As a
result, the most recent features have not yet been incorporated into the
version available on the CRAN repository. The current CRAN version
supports only the `linda` and `linda.plot` functions. For users who
require a broader range of functionalities, especially those related to
the analysis of longitudinal data, it is advisable to install the
development version directly from GitHub. This process necessitates the
prior installation of the `devtools` package.

``` r
install.packages("devtools")
```

Once `devtools` is installed, you can install `MicrobiomeStat` from
GitHub using the following command:

``` r
devtools::install_github("cafferychen777/MicrobiomeStat")
```

# Table of Contents

1.  [Citations](#citations)
    - [General Citation for MicrobiomeStat](#general-citation)
    - [Citation for Specialized `MicrobiomeStat`
      Functions](#specialized-citation)
    - [Important Note on CRAN Version](#cran-version-note)
2.  [Online Tutorials](#online-tutorials)
    - [Acquaint Yourself for a Seamless Experience](#acquaint-yourself)
    - [Discovering MicrobiomeStat](#explore-microbiomestat)
      - [Acknowledgements](#acknowledgements)
3.  [Benefits of Using MicrobiomeStat](#why-choose-microbiomestat)
    - [User Support](#user-support)
    - [Ongoing Development](#ongoing-development)
    - [Collaborative Development](#collaborative-development)
    - [Conclusion and Features](#conclusion)
    - [Demo Reports](#demo-reports)
4.  [Assistance & Contact Information](#support-contact)
5.  [Engage in Our Discord Community](#discord-community)
6.  [Share and Connect](#share-and-connect)

# Online Tutorials

`MicrobiomeStat` provides a comprehensive suite of tools for microbiome
data analysis, encompassing a variety of functions from data input to
visualization.

To acquaint users with `MicrobiomeStat`, we offer an extensive online
tutorial on GitBook. The tutorial covers the following areas:

- Installation and Configuration Instructions
  - These guidelines help ensure that your setup is correctly configured
    and optimized.
- Analysis Demonstrations Based on Real-World Scenarios
  - These demonstrations provide practical insights and skills.
- Code Examples for Practice
  - These examples allow users to familiarize themselves with
    `MicrobiomeStat` coding practices.
- Guides for Interpreting Results and Creating Visualizations
  - These guides help users understand and effectively present their
    data.
- Answers to Frequently Asked Questions
  - This section provides quick solutions to common questions.

## Acquaint Yourself for a Seamless Experience

For a seamless experience with `MicrobiomeStat`, make the most of these
enriching resources:

[**📘 Explore MicrobiomeStat
Tutorials**](https://www.microbiomestat.wiki)

## Discovering MicrobiomeStat

The realm of microbiome research is intricate and continually advancing.
The analytical tools selected can play a crucial role in navigating
through the research journey. In this scenario, `MicrobiomeStat` aims to
be a supportive companion.

### Acknowledgements

We stand on the shoulders of giants with `MicrobiomeStat`, and our
heartfelt gratitude goes out to the diligent and brilliant developers of
the dependencies that our package relies on. Their remarkable efforts
have not only made our work possible but have also significantly
elevated the standards of computational tools available to the
scientific community:

- Core Dependencies:
  - R (\>= 3.5.0), rlang, tibble
- Imported Packages:
  - ggplot2, matrixStats, lmerTest, foreach, modeest, vegan, dplyr,
    pheatmap, tidyr, ggh4x, ape, GUniFrac, scales, stringr, rmarkdown,
    knitr, pander, tinytex
- Suggested Packages:
  - ggrepel, parallel, ggprism, aplot, philentropy, forcats, yaml,
    biomformat, Biostrings

Furthermore, we extend our deepest appreciation and respect to the
trailblazers in the microbiome research community who have created and
maintained the following remarkable tools. Their pioneering work has
laid down paths through the complex landscape of microbiome data
analysis, and we are truly honored to walk alongside:

- `microbiomeutilities`, `phyloseq`, `microbiomemarker`,
  `MicrobiomeAnalyst`, `microbiomeeco`, `EasyAmplicon`, `STAMP`,
  `qiime2`, and `MicrobiotaProcess`

Their contributions inspire us to continue improving and expanding the
capabilities of `MicrobiomeStat`, and we sincerely hope our humble
addition proves to be a useful complement to the incredible array of
tools already available to researchers.

### User Support

`MicrobiomeStat` is designed with users in mind. Comprehensive
[documentation and tutorials](https://www.microbiomestat.wiki/) are
available to assist both novice and experienced researchers. Before
posting a question or issue, we encourage users to [check previous
questions and
issues](https://github.com/cafferychen777/MicrobiomeStat/issues?q=is%3Aissue+is%3Aclosed)
to see if the topic has already been addressed.

In case you have specific comments or questions about a particular
function’s documentation and you find that the RStudio’s search box
leads to a 404 error, you can directly access the function’s
documentation at
<https://cafferychen777.github.io/MicrobiomeStat/reference/index.html>.

If your question or issue has not been previously addressed, feel free
to [open a new issue on
GitHub](https://github.com/cafferychen777/MicrobiomeStat/issues). We are
here to help you navigate any challenges you may encounter.

### Ongoing Development

The `MicrobiomeStat` tool is under continuous development to incorporate
user feedback and to keep up with advancements in the field. We are
pleased to announce that the Shiny interface for `MicrobiomeStat` is now
officially available. This interface provides an interactive platform
for microbiome data analysis.

The Shiny interface will be maintained and updated along with the main
package. The Shiny application can be accessed directly at [this
link](https://cafferyyang.shinyapps.io/MicrobiomeStat_Shiny/). For users
who prefer a local setup or require more customization, the Shiny
application files and instructions are available on its dedicated
[GitHub
repository](https://github.com/cafferychen777/MicrobiomeStat-Shiny).

### Collaborative Development

`MicrobiomeStat` is an open-source tool, and we highly value
contributions from the community. If you have suggestions, improvements,
or feedback for future development directions and feature additions,
[pull requests](https://github.com/cafferychen777/MicrobiomeStat/pulls)
are welcomed, and you can also share your ideas in the [discussion
area](https://github.com/cafferychen777/MicrobiomeStat/discussions) of
our GitHub repository. Engage with other community members and help us
make `MicrobiomeStat` an even more useful tool for microbiome research.

### Conclusion

`MicrobiomeStat` aims to serve as a dependable and efficient resource
for microbiome data analysis. We extend an invitation to those who value
open-source collaboration to join our community and contribute to its
ongoing development.

| Feature                                                                                                                                                                                 | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [Data Import and Conversion](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object) | Accommodates multiple input formats from platforms such as [QIIME2](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object/importing-data-from-qiime2-into-microbiomestat), [Mothur](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object/fetching-data-from-mothur-into-microbiomestat), [DADA2](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object/integrating-data-from-dada2-into-microbiomestat), [Phyloseq](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object/navigating-data-from-phyloseq-into-microbiomestat) and others. |
| [Cross-sectional Study Analysis](https://www.microbiomestat.wiki/cross-sectional-study-design/unraveling-cross-sectional-studies-with-microbiomestat)                                   | Offers thorough analysis for cross-sectional studies.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| [Paired Sample Analysis](https://www.microbiomestat.wiki/paired-samples-analysis/unveiling-paired-samples-analysis-a-comprehensive-guide)                                               | Provides tools for paired samples analysis.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| [Longitudinal Study Analysis](https://www.microbiomestat.wiki/longitudinal-study-design/grasping-longitudinal-studies-introduction-and-dataset-overview)                                | Facilitates exploration of the temporal dynamics of the microbiome.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| Report Generation Functions                                                                                                                                                             | Includes individual report functions for [cross-sectional](https://www.microbiomestat.wiki/cross-sectional-study-design/cross-sectional-reporting-microbial-analysis-reports-with-microbiomestat), [paired](https://www.microbiomestat.wiki/paired-samples-analysis/automated-reporting-for-paired-studies-microbiomestats-integrated-analysis-reports), [longitudinal](https://www.microbiomestat.wiki/longitudinal-study-design/longitudinal-reporting-microbiome-analysis-automation-with-microbiomestat) study designs. Development of a Shiny interface for one-click reporting is underway.                                                                                                                                                                                                                                                                                                                                                         |
| Visualization Capabilities                                                                                                                                                              | Supports a broad range of visualization styles.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| Ongoing Development                                                                                                                                                                     | Committed to continuous refinement of existing features and addition of new functionalities.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |

This approach ensures that users can effortlessly navigate to the
specific sections of the `MicrobiomeStat` documentation, garnering
detailed information and guidelines for diverse analysis types. The
structure and accessibility assist users in leveraging `MicrobiomeStat`
effectively for their microbiome data analysis needs.

### Demo Reports

For those interested in seeing `MicrobiomeStat` in action, we have
prepared demo reports tailored to different study designs:

- [Cross-sectional Study Design: Reporting Microbial Analysis with
  MicrobiomeStat](https://www.microbiomestat.wiki/cross-sectional-study-design/cross-sectional-reporting-microbial-analysis-reports-with-microbiomestat)
- [Paired Samples Analysis: Reporting Microbial Analysis with
  MicrobiomeStat](https://www.microbiomestat.wiki/paired-samples-analysis/automated-reporting-for-paired-studies-microbiomestats-integrated-analysis-reports)
- [Longitudinal Study Design: Microbiome Analysis Automation with
  MicrobiomeStat](https://www.microbiomestat.wiki/longitudinal-study-design/longitudinal-reporting-microbiome-analysis-automation-with-microbiomestat)

We encourage you to explore these examples and discover the powerful
capabilities of our tool.

## Assistance & Contact Information

For assistance or inquiries, feel free to reach out to:

| Name                                                                         | Email                       |
|------------------------------------------------------------------------------|-----------------------------|
| [Dr. Jun Chen](https://scholar.google.com/citations?user=gonDvdwAAAAJ&hl=en) | <Chen.Jun2@mayo.edu>        |
| [Chen Yang](https://cafferyyang.owlstown.net/)                               | <cafferychen7850@gmail.com> |

## Engage in Our Discord Community

Join our Discord community to stay on top of the latest updates,
developments, and enhancements in `MicrobiomeStat`. Be part of vibrant
discussions, ask questions, share insights, and avail support from peers
and experts alike:

[Join the MicrobiomeStat Discord Server!](https://discord.gg/BfNvTJAt)

In our Discord server, an automated bot keeps you informed about every
package and tutorial update, ensuring you never miss out on new
features, improvements, and learning materials. Our active community
thrives on collaboration, feedback, and continuous learning, making it
an invaluable space for both novice and experienced researchers
navigating the world of microbiome data analysis. Stay connected, stay
informed, and let’s advance the field of microbiome data analysis
together!

## Share and Connect

Spread the word about `MicrobiomeStat` and stay connected through
various platforms!

[![Twitter](https://img.shields.io/twitter/url?url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&style=social)](https://twitter.com/intent/tweet?url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&text=Check%20out%20this%20awesome%20package%20for%20comprehensive%20and%20longitudinal%20microbiome%20analysis%20in%20R!)

[![Facebook](https://img.shields.io/badge/Share_on-Facebook-1877F2?logo=facebook&style=social)](https://www.facebook.com/sharer/sharer.php?u=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&quote=Check%20out%20this%20awesome%20package%20for%20comprehensive%20and%20longitudinal%20microbiome%20analysis%20in%20R!)

[![LinkedIn](https://img.shields.io/badge/Share_on-LinkedIn-0077B5?logo=linkedin&style=social)](https://www.linkedin.com/shareArticle?mini=true&url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&title=Check%20out%20this%20awesome%20package%20for%20comprehensive%20and%20longitudinal%20microbiome%20analysis%20in%20R!)

[![Reddit](https://img.shields.io/badge/Share_on-Reddit-FF4500?logo=reddit&style=social)](https://www.reddit.com/submit?url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&title=Check%20out%20this%20awesome%20package%20for%20comprehensive%20and%20longitudinal%20microbiome%20analysis%20in%20R!)

[![WhatsApp](https://img.shields.io/badge/Share_on-WhatsApp-25D366?logo=whatsapp&style=social)](https://wa.me/?text=Check%20out%20this%20awesome%20package%20for%20comprehensive%20and%20longitudinal%20microbiome%20analysis%20in%20R!%20https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat)

[![Slack](https://img.shields.io/badge/Share_on-Slack-4A154B?logo=slack&style=social)](https://slack.com/intl/en-cn/)

[![Email](https://img.shields.io/badge/Share_on-Gmail-D14836?logo=gmail&style=social)](mailto:?subject=Check%20out%20this%20awesome%20package%20for%20comprehensive%20and%20longitudinal%20microbiome%20analysis%20in%20R!&body=I%20found%20this%20amazing%20R%20package%20for%20microbiome%20analysis%20called%20%60MicrobiomeStat%60.%20Check%20it%20out%20here:%20https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat)
