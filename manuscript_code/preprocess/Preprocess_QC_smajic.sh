## The Smajic et al. data object was processed with scRNAbox: https://github.com/neurobioinfo/scrnabox

############################################################################ Step 0
############################################################################
############################################################################
############################################################################
screen -S scrnabx_cont

export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.03/scrnabox.pip
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning/scRNAbox_smajic/scRNAbox
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 0 --method SCRNA --container TRUE

############################################################################ Step 2
############################################################################
############################################################################
############################################################################

export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.03/scrnabox.pip
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning/scRNAbox_smajic/scRNAbox/
cd /home/fiorini9/scratch/machine_learning/scRNAbox_smajic/scRNAbox/
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 2


############################################################################ Step 3
############################################################################
############################################################################
############################################################################

export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.03/scrnabox.pip
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning/scRNAbox_smajic/scRNAbox/
cd /home/fiorini9/scratch/machine_learning/scRNAbox_smajic/scRNAbox/
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 3

############################################################################ Step 4
############################################################################
############################################################################
############################################################################

export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.03/scrnabox.pip
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning/scRNAbox_smajic/scRNAbox/
cd /home/fiorini9/scratch/machine_learning/scRNAbox_smajic/scRNAbox/
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 4


############################################################################ Step 5
############################################################################
############################################################################
############################################################################

export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.03/scrnabox.pip
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning/scRNAbox_smajic/scRNAbox/
cd /home/fiorini9/scratch/machine_learning/scRNAbox_smajic/scRNAbox/
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 5

############################################################################ Step 6
############################################################################
############################################################################
############################################################################

export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.03/scrnabox.pip
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning/scRNAbox_smajic/scRNAbox/
cd /home/fiorini9/scratch/machine_learning/scRNAbox_smajic/scRNAbox/
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 6

############################################################################ Step 7
############################################################################
############################################################################
############################################################################

export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.03/scrnabox.pip
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning/scRNAbox_smajic/scRNAbox/
cd /home/fiorini9/scratch/machine_learning/scRNAbox_smajic/scRNAbox/
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 7