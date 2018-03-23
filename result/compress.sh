# tar and compress every folder
find . -maxdepth 1 -mindepth 1 -type d -exec tar cvf {}.tar {} --remove-files  \;
pxz -v -0 ./*.tar

# remove empty folders
find . -maxdepth 1 -mindepth 1 -type d -empty -exec rmdir {}  \;
