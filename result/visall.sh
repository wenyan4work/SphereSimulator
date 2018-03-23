rm -rf ./ResultAll
mkdir ./ResultAll
find $PWD -maxdepth 2 -mindepth 2 -type f -name '*.vtp' -exec ln -s {} $PWD/ResultAll/ \;
find $PWD -maxdepth 2 -mindepth 2 -type f -name '*.pvtp' -exec ln -s {} $PWD/ResultAll/ \;
find $PWD -maxdepth 2 -mindepth 2 -type f -name '*.vtr' -exec ln -s {} $PWD/ResultAll/ \;
find $PWD -maxdepth 2 -mindepth 2 -type f -name '*.vtu' -exec ln -s {} $PWD/ResultAll/ \;
find $PWD -maxdepth 2 -mindepth 2 -type f -name '*.pvtu' -exec ln -s {} $PWD/ResultAll/ \;

# -exec tar cvf {}.tar {} \;
