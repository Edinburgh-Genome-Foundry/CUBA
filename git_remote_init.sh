# ./init_remote_git.sh root@139.59.162.184 CUBA prodserver
#
# Executes:
# ssh root@139.59.162.184 'git init CUBA.git/;cd CUBA.git;git config receive.denyCurrentBranch updateInstead'
# git remote add prodserver root@139.59.162.184:CUBA.git
echo "git init $2.git/"
ssh $1 "git init $2.git/;cd $2.git/;git config receive.denyCurrentBranch updateInstead"
git remote add $3 $1:$2.git
git push $3 master
