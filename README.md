## Archsphere
Point Clustering for geographic data

#For ArcMAP

Save ArchSphere.tbx and archsphereGIS.py on your local drive and add the tbx file to your ArcToolbox.
Double click the archsphere script to start the tool and set your parameters:

The input shapefile should be projected (in meter), best is UTM
You can get height infromation from a DEM using the Spatial Analyst tool "Extract Multi Values to Points" if you have not collected this information.
The Multiplier increments in 0.1 steps
The output location should be an existing folder.

1. Add your input shapefile
2. Set your height field. 
3. Set your weight field
4. Set the low Minimum Total Points
5. Set the high Minimum Total Points
6. Set the low Multiplier
7. Set the high Multiplier
8. Your output location  


