# Known issues
1. If one or more rank has no spheres, the VTP/VTU xml file output with an empty content will not be recognized by Paraview. 
2. If one or more rank has no spheres, the neighbor detection part may crash.
3. If one or more rank has no spheres, the initial partition() to broadcast and split spheres from root rank to all ranks runs fine. 