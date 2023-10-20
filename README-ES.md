MicrobiomeStat: Apoyo al análisis longitudinal del microbioma en R
================

<p align="center" style="margin:0; padding:0;">
<img src="man/figures/logo.png" alt="Logo de MicrobiomeStat" width="400" style="margin:0; padding:0;"/>
</p>
<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/MicrobiomeStat)](https://cran.r-project.org/package=MicrobiomeStat)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/MicrobiomeStat)](https://cran.r-project.org/package=MicrobiomeStat)
[![Licencia: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

[EN](README.md) \| [CN](README-CN.md) \| [ES](README-ES.md)

El paquete `MicrobiomeStat` es una herramienta de R dedicada para
explorar **datos longitudinales del microbioma**. También acomoda datos
multi-ómicos y estudios transversales, valorando los esfuerzos
colectivos dentro de la comunidad. Esta herramienta tiene como objetivo
apoyar a los investigadores a través de sus extensas consultas
biológicas a lo largo del tiempo, con un espíritu de gratitud hacia los
recursos existentes de la comunidad y un ethos de colaboración para
avanzar en la investigación del microbioma.

# Noticias

📢 **Actualización** (20 de octubre): La interfaz de Shiny ya está
oficialmente disponible para su uso. Actualmente está configurada para
manejar análisis de conjuntos de datos pequeños a medianos. Puede
acceder a la interfaz a través de [este
enlace](https://a95dps-caffery-chen.shinyapps.io/MicrobiomeStat_Shiny/).

En caso de que las limitaciones del servidor afecten su análisis, o para
aquellos que prefieren trabajar con módulos más pequeños, recomendamos
usar nuestro paquete directamente. Para una implementación más flexible,
considere clonar nuestro repositorio Shiny desde
[aquí](https://github.com/cafferychen777/MicrobiomeStat-Shiny) y
desplegarlo en su servidor local o computadora.

Agradecemos su comprensión y participación continua.

# Citas

## Cita general para MicrobiomeStat

Si está utilizando funciones más allá de `linda` y `linda.plot`, cite de
la siguiente manera, hasta que se publique una versión preimpresa:

    @Manual{,
      title = {MicrobiomeStat: Métodos Estadísticos y de Visualización Completos para Datos de Microbioma y Multi-Ómicos},
      author = {Xianyang Zhang y Jun Chen y Caffery (Chen) Yang},
      year = {2023},
      note = {Versión del paquete R 1.1.1},
      url = {https://www.microbiomestat.wiki},
    }

## Cita para funciones especializadas de `MicrobiomeStat`

Si está utilizando las funciones `linda`, `linda.plot`,
`generate_taxa_association_test_long`, `generate_taxa_test_pair`,
`generate_taxa_test_single`, o `generate_taxa_trend_test_long`, cite el
siguiente artículo:

    @article{zhou2022linda,
      title={LinDA: modelos lineales para el análisis de abundancia diferencial de datos composicionales del microbioma},
      author={Zhou, Huijuan y He, Kejun y Chen, Jun y Zhang, Xianyang},
      journal={Genome biology},
      volume={23},
      number={1},
      pages={1--23},
      year={2022},
      publisher={BioMed Central}
    }

Actualizaremos las pautas de citación tan pronto como se publique el
preprint.

## Nota importante sobre la versión CRAN

El paquete `MicrobiomeStat` está en desarrollo continuo. Como resultado,
las características más recientes aún no se han incorporado en la
versión disponible en el repositorio CRAN. La versión actual de CRAN
solo admite las funciones `linda` y `linda.plot`. Para los usuarios que
requieren un rango más amplio de funcionalidades, especialmente aquellas
relacionadas con el análisis de datos longitudinales, se recomienda
instalar la versión de desarrollo directamente desde GitHub. Este
proceso requiere la instalación previa del paquete `devtools`.

``` r
install.packages("devtools")
```

Una vez instalado `devtools`, puede instalar `MicrobiomeStat` desde
GitHub utilizando el siguiente comando:

``` r
devtools::install_github("cafferychen777/MicrobiomeStat")
```

# Tabla de contenidos

1.  [Citas](#citations)
    - [Cita general para MicrobiomeStat](#general-citation)
    - [Cita para funciones especializadas de
      `MicrobiomeStat`](#specialized-citation)
    - [Nota importante sobre la versión CRAN](#cran-version-note)
2.  [Tutoriales en línea](#online-tutorials)
    - [Familiarízate para una experiencia sin
      problemas](#acquaint-yourself)
    - [Descubriendo MicrobiomeStat](#explore-microbiomestat)
      - [Exploración de características](#feature-exploration)
      - [Agradecimientos](#acknowledgements)
3.  [Beneficios de usar MicrobiomeStat](#why-choose-microbiomestat)
    - [Soporte al usuario](#user-support)
    - [Desarrollo en curso](#ongoing-development)
    - [Desarrollo colaborativo](#collaborative-development)
    - [Conclusión y características](#conclusion)
    - [Informes de demostración](#demo-reports)
4.  [Asistencia e información de contacto](#support-contact)
5.  [Participa en nuestra comunidad de Discord](#discord-community)
6.  [Compartir y conectar](#share-and-connect)

# Tutoriales en línea

`MicrobiomeStat` proporciona un conjunto completo de herramientas para
el análisis de datos del microbioma, que abarca una variedad de
funciones desde la entrada de datos hasta la visualización.

Para familiarizar a los usuarios con `MicrobiomeStat`, ofrecemos un
extenso tutorial en línea en GitBook. El tutorial cubre las siguientes
áreas:

- Instrucciones de instalación y configuración
  - Estas pautas ayudan a garantizar que su configuración esté
    correctamente configurada y optimizada.
- Demostraciones de análisis basadas en escenarios del mundo real
  - Estas demostraciones proporcionan conocimientos prácticos y
    habilidades.
- Ejemplos de código para la práctica
  - Estos ejemplos permiten a los usuarios familiarizarse con las
    prácticas de codificación de `MicrobiomeStat`.
- Guías para interpretar resultados y crear visualizaciones
  - Estas guías ayudan a los usuarios a entender y presentar eficazmente
    sus datos.
- Respuestas a preguntas frecuentes
  - Esta sección proporciona soluciones rápidas a preguntas comunes.

## Familiarízate para una experiencia sin problemas

Para una experiencia sin problemas con `MicrobiomeStat`, aprovecha al
máximo estos recursos enriquecedores:

[**📘 Explora los tutoriales de
MicrobiomeStat**](https://www.microbiomestat.wiki)

## Descubriendo MicrobiomeStat

El ámbito de la investigación del microbioma es intrincado y avanza
continuamente. Las herramientas analíticas seleccionadas pueden
desempeñar un papel crucial en la navegación a través del viaje de
investigación. En este escenario, `MicrobiomeStat` tiene como objetivo
ser un compañero de apoyo.

### Exploración de características

Para entender mejor las capacidades de `MicrobiomeStat`, hemos esbozado
sus características y funcionalidades en nuestro sitio web, junto con
algunas referencias contextuales a otras herramientas para una
perspectiva más informada:

- [Explorando características de análisis
  longitudinal](https://www.microbiomestat.wiki/introduction/unveiling-microbiomestat-a-glimpse-into-diverse-microbiome-analysis-solutions/exploring-longitudinal-analysis-features-microbiomestat-among-other-packages)

- [Explorando características de análisis
  integrado](https://www.microbiomestat.wiki/introduction/unveiling-microbiomestat-a-glimpse-into-diverse-microbiome-analysis-solutions/exploring-feature-offerings-microbiomestat-among-integrated-analysis-packages)

### Agradecimientos

Nos apoyamos en los hombros de gigantes con `MicrobiomeStat`, y nuestra
gratitud sincera va a los diligentes y brillantes desarrolladores de las
dependencias en las que se basa nuestro paquete. Sus esfuerzos notables
no solo han hecho posible nuestro trabajo, sino que también han elevado
significativamente los estándares de las herramientas computacionales
disponibles para la comunidad científica:

- Dependencias principales:
  - R (\>= 3.5.0), rlang, tibble
- Paquetes importados:
  - ggplot2, matrixStats, lmerTest, foreach, modeest, vegan, dplyr,
    pheatmap, tidyr, ggh4x, ape, GUniFrac, scales, stringr, rmarkdown,
    knitr, pander, tinytex
- Paquetes sugeridos:
  - ggrepel, parallel, ggprism, aplot, philentropy, forcats, yaml,
    biomformat, Biostrings

Además, extendemos nuestro más profundo agradecimiento y respeto a los
pioneros en la comunidad de investigación del microbioma que han creado
y mantenido las siguientes herramientas notables. Su trabajo pionero ha
trazado caminos a través del complejo paisaje del análisis de datos del
microbioma, y nos sentimos verdaderamente honrados de caminar junto a
ellos:

- `microbiomeutilities`, `phyloseq`, `microbiomemarker`,
  `MicrobiomeAnalyst`, `microbiomeeco`, `EasyAmplicon`, `STAMP`,
  `qiime2`, y `MicrobiotaProcess`

Sus contribuciones nos inspiran a continuar mejorando y expandiendo las
capacidades de `MicrobiomeStat`, y esperamos sinceramente que nuestra
humilde adición resulte ser un complemento útil para la increíble
variedad de herramientas ya disponibles para los investigadores.

### Soporte al usuario

`MicrobiomeStat` está diseñado pensando en los usuarios. Hay
[documentación y tutoriales](https://www.microbiomestat.wiki/)
disponibles para ayudar tanto a los investigadores novatos como a los
experimentados. Antes de publicar una pregunta o problema, alentamos a
los usuarios a [verificar las preguntas y problemas
anteriores](https://github.com/cafferychen777/MicrobiomeStat/issues?q=is%3Aissue+is%3Aclosed)
para ver si el tema ya ha sido abordado.

En caso de que tenga comentarios o preguntas específicas sobre la
documentación de una función en particular y encuentre que el cuadro de
búsqueda de RStudio conduce a un error 404, puede acceder directamente a
la documentación de la función en
<https://cafferychen777.github.io/MicrobiomeStat/reference/index.html>.

Si su pregunta o problema no ha sido abordado previamente, no dude en
[abrir un nuevo problema en
GitHub](https://github.com/cafferychen777/MicrobiomeStat/issues).
Fomentamos activamente el intercambio y la colaboración entre
investigadores, desarrolladores y usuarios de todo el mundo. Para
promover una comunicación más fluida y eficiente a nivel global, le
sugerimos hacer sus preguntas en inglés. Esto permite que un público más
amplio comprenda y participe en la discusión, y nos asegura poder
brindarle asistencia de manera más rápida y precisa. ¡Gracias por su
comprensión y apoyo! Estamos aquí para ayudarlo a navegar cualquier
desafío que pueda encontrar.

### Desarrollo en curso

Asegurar que `MicrobiomeStat` siga siendo una herramienta líder en su
categoría requiere un desarrollo continuo. Estamos dedicados a las
actualizaciones regulares y a abordar los comentarios de los usuarios.
Como parte de nuestros esfuerzos de mejora continua, hemos desarrollado
una interfaz Shiny para `MicrobiomeStat`, que ofrece una forma
interactiva y fácil de usar para realizar análisis de datos de
microbiomas. La interfaz Shiny se mantiene y mejora activamente junto
con el paquete principal. Puedes acceder a los archivos de la aplicación
Shiny y las instrucciones para la configuración local en su [repositorio
de GitHub](https://github.com/cafferychen777/MicrobiomeStat-Shiny)
dedicado.

### Desarrollo colaborativo

`MicrobiomeStat` es una herramienta de código abierto, y valoramos mucho
las contribuciones de la comunidad. Si tiene sugerencias, mejoras o
comentarios para futuras direcciones de desarrollo y adiciones de
características, se aceptan [solicitudes de
extracción](https://github.com/cafferychen777/MicrobiomeStat/pulls), y
también puede compartir sus ideas en el [área de
discusión](https://github.com/cafferychen777/MicrobiomeStat/discussions)
de nuestro repositorio de GitHub. Participa con otros miembros de la
comunidad y ayúdanos a hacer de `MicrobiomeStat` una herramienta aún más
útil para la investigación del microbioma.

### Conclusión

`MicrobiomeStat` tiene como objetivo servir como un recurso confiable y
eficiente para el análisis de datos del microbioma. Extendemos una
invitación a aquellos que valoran la colaboración de código abierto para
unirse a nuestra comunidad y contribuir a su desarrollo continuo.

| Característica                                                                                                                                                                                 | Descripción                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [Importación y conversión de datos](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object) | Acomoda múltiples formatos de entrada de plataformas como [QIIME2](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object/importing-data-from-qiime2-into-microbiomestat), [Mothur](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object/fetching-data-from-mothur-into-microbiomestat), [DADA2](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object/integrating-data-from-dada2-into-microbiomestat), [Phyloseq](https://www.microbiomestat.wiki/setting-up-microbiomestat-installation-and-data-preparation/laying-the-foundation-creating-the-microbiomestat-data-object/navigating-data-from-phyloseq-into-microbiomestat) y otros. |
| [Análisis de estudios transversales](https://www.microbiomestat.wiki/cross-sectional-study-design/unraveling-cross-sectional-studies-with-microbiomestat)                                      | Ofrece un análisis exhaustivo para estudios transversales.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| [Análisis de muestras emparejadas](https://www.microbiomestat.wiki/paired-samples-analysis/unveiling-paired-samples-analysis-a-comprehensive-guide)                                            | Proporciona herramientas para el análisis de muestras emparejadas.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| [Análisis de estudios longitudinales](https://www.microbiomestat.wiki/longitudinal-study-design/grasping-longitudinal-studies-introduction-and-dataset-overview)                               | Facilita la exploración de la dinámica temporal del microbioma.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| Funciones de generación de informes                                                                                                                                                            | Incluye funciones de informes individuales para [estudios transversales](https://www.microbiomestat.wiki/cross-sectional-study-design/cross-sectional-reporting-microbial-analysis-reports-with-microbiomestat), [emparejados](https://www.microbiomestat.wiki/paired-samples-analysis/automated-reporting-for-paired-studies-microbiomestats-integrated-analysis-reports), [longitudinales](https://www.microbiomestat.wiki/longitudinal-study-design/longitudinal-reporting-microbiome-analysis-automation-with-microbiomestat) diseños de estudio. Se está desarrollando una interfaz Shiny para informes con un solo clic.                                                                                                                                                                                                                                                                                                                        |
| Capacidades de visualización                                                                                                                                                                   | Admite una amplia gama de estilos de visualización.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Desarrollo en curso                                                                                                                                                                            | Comprometido con el refinamiento continuo de las características existentes y la adición de nuevas funcionalidades.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |

Este enfoque asegura que los usuarios puedan navegar sin esfuerzo por
las secciones específicas de la documentación de `MicrobiomeStat`,
obteniendo información detallada y directrices para diversos tipos de
análisis. La estructura y accesibilidad ayudan a los usuarios a
aprovechar eficazmente `MicrobiomeStat` para sus necesidades de análisis
de datos del microbioma.

### Informes de demostración

Para aquellos interesados en ver `MicrobiomeStat` en acción, hemos
preparado informes de demostración adaptados a diferentes diseños de
estudio:

- [Diseño de estudio transversal: Informes de análisis microbiano con
  MicrobiomeStat](https://www.microbiomestat.wiki/cross-sectional-study-design/cross-sectional-reporting-microbial-analysis-reports-with-microbiomestat)
- [Análisis de muestras emparejadas: Informes de análisis microbiano con
  MicrobiomeStat](https://www.microbiomestat.wiki/paired-samples-analysis/automated-reporting-for-paired-studies-microbiomestats-integrated-analysis-reports)
- [Diseño de estudio longitudinal: Automatización del análisis del
  microbioma con
  MicrobiomeStat](https://www.microbiomestat.wiki/longitudinal-study-design/longitudinal-reporting-microbiome-analysis-automation-with-microbiomestat)

Te animamos a explorar estos ejemplos y descubrir las poderosas
capacidades de nuestra herramienta.

## Asistencia e información de contacto

Para asistencia o consultas, no dudes en contactar a:

| Nombre                                                                       | Correo electrónico          |
|------------------------------------------------------------------------------|-----------------------------|
| [Dr. Jun Chen](https://scholar.google.com/citations?user=gonDvdwAAAAJ&hl=en) | <Chen.Jun2@mayo.edu>        |
| [Chen Yang](https://cafferyyang.owlstown.net/)                               | <cafferychen7850@gmail.com> |

## Participa en nuestra comunidad de Discord

Únete a nuestra comunidad de Discord para mantenerte al día con las
últimas actualizaciones, desarrollos y mejoras en `MicrobiomeStat`. Sé
parte de vibrantes discusiones, haz preguntas, comparte ideas y obtén
apoyo de compañeros y expertos por igual:

[¡Únete al servidor de Discord de
MicrobiomeStat!](https://discord.gg/BfNvTJAt)

En nuestro servidor de Discord, un bot automatizado te mantiene
informado sobre cada actualización del paquete y tutorial, asegurándote
de que nunca te pierdas las nuevas características, mejoras y materiales
de aprendizaje. Nuestra comunidad activa prospera en la colaboración, la
retroalimentación y el aprendizaje continuo, convirtiéndola en un
espacio invaluable tanto para los investigadores novatos como para los
experimentados que navegan por el mundo del análisis de datos del
microbioma. ¡Mantente conectado, mantente informado y avancemos juntos
en el campo del análisis de datos del microbioma!

## Compartir y conectar

¡Difunde la palabra sobre `MicrobiomeStat` y mantente conectado a través
de varias plataformas!

[![Twitter](https://img.shields.io/twitter/url?url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&style=social)](https://twitter.com/intent/tweet?url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&text=¡Echa%20un%20vistazo%20a%20este%20increíble%20paquete%20para%20el%20análisis%20completo%20y%20longitudinal%20del%20microbioma%20en%20R!)

[![Facebook](https://img.shields.io/badge/Compartir_en-Facebook-1877F2?logo=facebook&style=social)](https://www.facebook.com/sharer/sharer.php?u=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&quote=¡Echa%20un%20vistazo%20a%20este%20increíble%20paquete%20para%20el%20análisis%20completo%20y%20longitudinal%20del%20microbioma%20en%20R!)

[![LinkedIn](https://img.shields.io/badge/Compartir_en-LinkedIn-0077B5?logo=linkedin&style=social)](https://www.linkedin.com/shareArticle?mini=true&url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&title=¡Echa%20un%20vistazo%20a%20este%20increíble%20paquete%20para%20el%20análisis%20completo%20y%20longitudinal%20del%20microbioma%20en%20R!)

[![Reddit](https://img.shields.io/badge/Compartir_en-Reddit-FF4500?logo=reddit&style=social)](https://www.reddit.com/submit?url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat&title=¡Echa%20un%20vistazo%20a%20este%20increíble%20paquete%20para%20el%20análisis%20completo%20y%20longitudinal%20del%20microbioma%20en%20R!)

[![WhatsApp](https://img.shields.io/badge/Compartir_en-WhatsApp-25D366?logo=whatsapp&style=social)](https://wa.me/?text=¡Echa%20un%20vistazo%20a%20este%20increíble%20paquete%20para%20el%20análisis%20completo%20y%20longitudinal%20del%20microbioma%20en%20R!%20https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat)

[![Slack](https://img.shields.io/badge/Compartir_en-Slack-4A154B?logo=slack&style=social)](https://slack.com/intl/en-cn/)

[![Email](https://img.shields.io/badge/Share_on-Gmail-D14836?logo=gmail&style=social)](mailto:?subject=Check%20out%20this%20awesome%20package%20for%20comprehensive%20and%20longitudinal%20microbiome%20analysis%20in%20R!&body=I%20found%20this%20amazing%20R%20package%20for%20microbiome%20analysis%20called%20%60MicrobiomeStat%60.%20Check%20it%20out%20here:%20https%3A%2F%2Fgithub.com%2Fcafferychen777%2FMicrobiomeStat)
