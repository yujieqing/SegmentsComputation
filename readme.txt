This code is for segments computation during the process of voxel-based ionosphere tomography inversion. Both ray tracing and tradition method of segments compuation were implemented. The code is tested under VS2017 compiler. 
The output segments can be visualized by paraview or VisIT software after conversion using scripts 'toVtp.py'. To use the script, python3.7 is required and vtk libary should be installed.  
./
	/main.cpp: the main function
	/Common.h,Common.cpp: common functions for both methods
	/tradition.h,tradition.cpp: tradition method 
	/raytracing.h,raytracing.cpp: ray tracing method
	/pugixml.hpp,pugixml.cpp,pugiconfig.hpp: a xml file libary. These files are downloaded from https://pugixml.org.
	/toVtp.py, vtkXMLFile.py: scripts for converting the resultant segments into .vtp format, which can be recognized by paraview(https://www.paraview.org/) or VisIT (https://wci.llnl.gov/simulation/computer-codes/visit/) software.
	
	/input
		/rayCoordinates.txt：sample data for input ray coordinates
		/rayCoordinates10.txt: a sample data with 10 rays taken from rayCoordinates.txt
		
		/voxelModels 
			ten different voxel models
	
	/output
		time_cost.txt: time cost for both two methods
		tradition_sparse(i).txt: the output of tradition method in sparse matrix format
		raytracing_sparse(i).txt: the output of raytracing method in sparse matrix format
		raytracing_centSize(i).txt: the output of raytracing method in centSize format, which can be read and converted into .vtp by toVtp.py
		tradition_centSize(i).txt: the output of tradition method in centSize format, which can be read and converted into .vtp by toVtp.py
		
For any question, please contact Dr. Jieqing YU (yujieqing@cumt.edu.cn)