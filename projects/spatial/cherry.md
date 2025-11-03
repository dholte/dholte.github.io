---
layout: default
title: Hunting Cherries
---

# Hunting Cherries

## Background 

NOTE: Typically, spatial analysis consists of five key stages: understanding your goal, preparing data, choosing suitable tools and techniques, performing the research, and estimating results.

## Analysis Goal

Identify locations within the City of Edmonton most likely to contain large quantities of unpicked, ripe saskatoon berries.

## Initial Data Preparation

The first dataset used will be the City of Edmonton Trees dataset which includes approximately 470,000 trees. Thankfully, a subset of this data including only the 59,392 edible fruit trees is also available on the open data portal. The data includes the following 19 columns of attributes for each tree:

*IMAGE PLACEHOLDER datacolumns*

A master copy of the unaltered data will be saved, and a working copy will be used. The following columns are clearly not useful for this analysis and will be removed from the working copy: Bears Edible Fruit (non fruit bearing trees are already excluded from this dataset), COUNT (irrelevant), and OWNER (all trees in this dataset are owned by City of Edmonton Parks).

Next, SPECIES_BOTANICAL, GENUS, and SPECIES can all be removed, so long as all or most of the trees have data in the SPECIES_COMMON. 

*IMAGE PLACEHOLDER dataprep1*

Their removal is validated after using the COUNTA function to confirm 100% of trees have their common species name included. The Point Location column contains the location of each tree in Well-Known Text (WKT) format, so the LATITUDE, LONGITUDE, AND LOCATION columns can all be removed. 

To improve performance and reduce processing time, the data will be filtered to include only cherry trees before it is imported into QGIS. 

## Query Builder 
 





To begin the analysis, it is assumed that the following characteristics will *increase* the likelihood of a specific location having large quantities of unpicked saskatoon berries:

1. High concentration of saskatoon trees. No justification required.
2. High concentration of healthy saskatoon trees. Healthy trees bear larger amounts of high quality fruit in comparison to non-healthy trees.
3. High concentration of saskatoon trees planted 6-8 years ago. Although saskatoons begin producing fruit after 3-5 years, ideal production occurs between 6-8 years of age (Source: University of Saskatchewan College of Agriculture and Bioresources).


NOTE: saskatoon-juniper rust, keep at least 2 km away from native juniper trees to avoid. productivity can be limited by shade
