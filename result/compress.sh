# tar and compress every folder
# find . -maxdepth 1 -mindepth 1 -type d -exec tar cvf {}.tar {} --remove-files  \;
# pxz -v -0 ./*.tar

# 7z tar and compress
find . -maxdepth 1 -mindepth 1 -type d -exec 7z a -mx=5 -mmt=on -sdel {}.7z {} \;

# remove empty folders
find . -maxdepth 1 -mindepth 1 -type d -empty -exec rmdir {}  \;
