#docker pull memesuite/memesuite:latest

instanceDir=/Users/javrodher/Work/RStudio-PRJs/TCGA_PanCan/data/memeInstance
mkdir $instanceDir
mkdir ${instanceDir}/tmp

cp /Users/javrodher/Work/biodata/genomes/hg38/hg38.fa ${instanceDir}/tmp/
cp /Users/javrodher/Work/biodata/JASPAR/homo_sapiens/CTCF_MA0139.1.meme ${instanceDir}/tmp/

docker run -d -v `pwd`:/home/meme -v ${instanceDir}/tmp:/tmp --name meme-container memesuite/memesuite:latest

# NOW GO TO DOCKER DESKTOP AND OPEN THE CLI terminal for the meme-container #
# Then RUN:

fimo /tmp/CTCF_MA0139.1.meme /tmp/hg38.fa
