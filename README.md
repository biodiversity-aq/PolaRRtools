# PolaRRtools

PolaRRtools is an R package develloped by the POLA3R-team (part of the LifeWatch Belgium infrastructure) to ease the re-use and analysis of molecular microbial datasets. This package focuses on finding and downloading sequence data from International Nucleotide Sequence Database Collaboration (INSDC) together with any associated metadata or environmental measurement data. PolaRRtools also contains tools to automize the quality controll and re-formatting in standardized formats of metadata, as well as other hand data manipulation tools.

## Introduction
Ecological research of microbes increasingly relies high throughput DNA sequencing techniques to answer questions on biodiversity, functional composition and biogeography. Common techniques include genome sequencing of isolated or cultured microorganisms, community-wide (i.e. 'meta-') sequencing of targeted marker genes (amplicon sequencing), up to sequencing all fragmented DNA or RNA in environmental samples to post-hoc re-assemple metagenomes or metatranscriptomes of the microbial communities (shotgun sequencig). The genomics or transcriptomics (in short: 'omics) datasets that are generated in this process are highly rich sources of information that stay relevant to the scientific community for years, perhaps even centuries, after the the results have been published. These combining these 'historical' molecular datasets and re-analyzing them can, for instance help to benchmark community composition changes agains past conditions, or provide a wider geographical or temporal scope for studies that are limited in resources. 

However, while it has become a common practice to submit molecular data (preferable annotated with metadata and associated environmental measurements) to online repositories (i.e. INSDC), it can be a remarkably complex venture to find and download the these data for secondary users. With the PolaRRtools package, we aim to make these data much more accessible. More specifically, PolaRRtools provides tools to download nucleotide sequence data from the INSDC databases directly through the R-console. Data can be downloaded in bulk (i.e. multiple projects or samples through a single command), and the functions can be piped to existing data processing pipelines, significanty reducting the time investment cost of data re-use. PolaRRtools also takes care of any associated metadata or environmental measurement data, allowing it to be downloaded, quality controlled and re-formatted in standardized formats, ready for use in R.

## Background
PolaRRtools R-package is part of the LifeWatch (Belgium) Virtual Research Environment (VRE). LifeWatch was established as part of the European Strategy Forum on Research Infrastructure (ESFRI) and aims to provide an infrasctructure that supports and facilitate biodiversity research, which is done by develloping products, like data portals, R-packages, and other tools for the scientific comunity.

As part of this framework, PolaRRtools was also made to closely integrate with the Polar’Omics Links in Arctic-Antarctic-Alpine (3 A's) Research portal, or POLA3R for short. POLA3R is a thematic online portal (also part of LifeWatch) that focusses on molecular biodiversity data resources from polar regions. On POLA3R molecular (public) datasets are enriched with the associated metadata and environmental information (often dislocated from the data in publication tables, apendices of a publication, or even local hard drives). On the POLA3R poral, we store the metadata and environmetal data, and link this to the associated publications and the sequences on INSDC. To allow interoperability with other systems, the portal is also designed to operate between different data archiving standards, such as the Minimum Information on any (x) Sequence (MIxS) as well as DarwinCore. Datasets that are listed on POLA3R are also registered on GBIF to increase their discoverability. As such, POLA3R aims to provide a hub for the polar scientific community, where they can discover high quality and complete molecular biodiversity data.

With PolaRRtools, we focues on the most used databases for sequence data that are part of the International Nucleotide Sequence Database Collaboration (INSDC). More specifically, the functions in PolaRRtools are made to retrieve data from the Sequence Read Archive (SRA) database of the National (American) Center for Biotechnolofy Information (NCBI), and the European Nucleotide Archive (ENA) database of the European Bioinformatics Institute (EBI, part of EMBL). These are stable, interlinked and goverment funded institutions that provide a long-term and safe solution to nucleotide data storage. 

## Additional Information
For more information, also see 
http://www.lifewatch.be
http://www.biodiversity.aq

## Metadata and environmental measurement data download
In the current version (2019) only focusses on data (sequence and metadata) from the International Nucleotide Sequence Database Collaboration (INSDC). The aim is to expand this offer in later updates.

To understand how to find and retrieve nucleotide sequence data, it is good to have a basic idea on the structure of data on  INSDC. On these databases, samples (also called 'Runs', the smalest unit of nucleotide sequence grouping) are grouped by biological origin into 'BioSamples', which in term is grouped into a project, called a BioProject. Each sample/run, BioSample and BioProject is designated a unique identifer number with a fixed prefix (that is: SRR for Run, SAMN for BioSample and PRJ for BioProject). Because the American, European and Japanese databases are interlinked, the same data can be downloaded through SRA as well as ENA.

In the PolaRRtools packages, we access data through the highest organisation level: the BioProject. This usually corresponds to a datasets that was analyzed in a (single) publication, thus combining all biological nucleotide sequence data related to a single initiative, originating from a single organization. 

PolaRRtools has two functions to access the metadata and any associated environmental measurements that were committed to INSDC alongside the sequence data. This accessory information can include a publication reference, geographical coordinates of the sample, sequence and lab protocol references, date and time of sampling, or other measurements such as conductivity, depth, pH, phosphate concentration,...

The get.BioProject.metadata.INSDC function will fetch the most basic metadata of a BioProject from the INSDC repositories using E-utils API function of NCBI. This does not require user identification to access NCBI with an API-key, but is very limited in the amount of metadata that can be retrieved. These basic metadata typically include Run number, release date, load date, spots, bases, av_MB and download path, but not any environmental measurements.

To access the complete (that is, all the information submitted to INSDC by the author) dataset of meta- and ancillory data, the get.sample.attributes.INSDC function will need to be used. This will fetch the all the metadata of the samples given to the sampleID or BioPrjct arguments, but requires a user-provided API-key to use the Entrez Programming Utilities (E-utilities) from NCBI.

See https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ for more info on why NCBI requires users to identify themselves with an API-key, and how to get one.


## Sequencedata download

Sequences can be downloaded using the download.sequences.INSDC function. This will collect all the nucleotide sequence data associated with the user-provided BioProject numbers, and will write those data to a user-provided path (computer location) as \*.fastq.gz files.

download.sequences.INSDC can also call get.BioProject.metadata.INSDC to simultaniously collect all the metadata associated with the sequence data.



