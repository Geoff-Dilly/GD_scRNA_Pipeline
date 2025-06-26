#!/bin/bash

cd "__HOME_DIR__";

R --quiet --no-restore --file="Scripts/Load_10x_Data.R";
#!R --quiet --no-restore --file="Scripts/Run_SoupX.R";
R --quiet --no-restore --file="Scripts/Run_DoubletFinder.R";
R --quiet --no-restore --file="Scripts/Merge_and_Normalize.R";
R --quiet --no-restore --file="Scripts/Cluster_and_ID_Cells.R"; 
R --quiet --no-restore --file="Scripts/Identify_Marker_Genes.R";
#!R --quiet --no-restore --file="Scripts/Rename_Clusters.R";
R --quiet --no-restore --file="Scripts/Make_Plots.R";
#!R --quiet --no-restore --file="Scripts/One_Var_DGE_Analysis.R";
