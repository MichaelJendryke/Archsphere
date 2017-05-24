## Archsphere
Point clustering for archeologically data

### Reference
Caspari, Gino; Jendryke, Michael (2017): Archsphere – A Cluster Algorithm for Archaeological Applications. Journal of Archaeological Science: Reports. [in review]

### How to use in ESRI/ArcMAP

Save ArchSphere.tbx and archsphereGIS.py to your local drive and add the tbx file to your ArcToolbox.
Double click the archsphere script to start the tool and set your parameters:

1. Add your input shapefile
2. Set your height field 
3. Set your weight field
4. Set the low Minimum Total Number of Points
5. Set the high Minimum Total Number of Points
6. Set the low Multiplier
7. Set the high Multiplier
8. Your output folder location  

#### Consider
- The input shapefile should be projected (in meter), best is UTM.
- Use Spatial Analyst tool "Extract Multi Values to Points" to get height information from a DEM or simply set the field to whatever you want if you have not collected height information.
- By setting different low and high values for Minimum Total Number of Points and Multiplier a loop is created and multiple shapes are created.  
- The Multiplier increments in 0.1 steps.
- The output location should be an existing folder.

### How to use the algorithm directly in python  
coming soon



