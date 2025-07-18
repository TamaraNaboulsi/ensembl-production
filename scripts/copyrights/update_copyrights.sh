#/bin/sh
# NOTE THIS script need `gh` command line tool available on https://github.com/cli/cli#installation
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2025] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License

usage () {
  echo "Usage: $0 [path_repo_list_file or repo full path] [reviewer] [assignee]"
}

if [ -z "$1" ] ; then
    usage
    exit 1
fi

if [ -f "$1" ]; then
  repositories=`cat ${1}`
else
  repositories=$1
fi

hash gh 2>/dev/null || {
  echo >&2 "This script required 'gh' library. See https://github.com/cli/cli#installation and rerun."
  exit 1
}

tmp_dir="$HOME/tmp"
mkdir -p $tmp_dir
year=`date +'%Y'`

for repo in $repositories; do
  # checkout default branch
  echo "--------------------"
  echo $repo
  if [[ $repo != Ensembl* ]]; then
    repo="Ensembl/$repo"
  fi
  rm -rf ${tmp_dir}/${repo}
  git clone --depth 1 --branch main git@github.com:${repo} ${tmp_dir}/${repo}
  if [ $? -eq 0 ]; then
    cd ${tmp_dir}/${repo}
    git push origin --delete bau/copyright-${year}
    git checkout -b bau/copyright-${year}
    perl ${ENSEMBL_ROOT_DIR}/ensembl/misc-scripts/annual_copyright_updater.sh
    git commit -a -m "${year} copyright update"
    if [ $? -eq 0 ]; then
      git push --set-upstream origin bau/copyright-${year}
      if [ $? -eq 0 ]; then
        if [ "$#" -eq 2 ]; then
          gh pr create --title "Annual copyright update ${year}" --body "${year} annual copyright file updates" --head bau/copyright-${year} --reviewer $2 --base main
        elif [ "$#" -eq 3 ]; then
          gh pr create --title "Annual copyright update ${year}" --body "${year} annual copyright file updates" --head bau/copyright-${year} --reviewer $2 --assignee $3 --base main
        else
          gh pr create --title "Annual copyright update ${year}" --body "${year} annual copyright file updates" --head bau/copyright-${year} --base main
        fi
      else
        echo 'failed to push commits and open a pull request.';
        git push origin --delete bau/copyright-${year}
      fi
    else
      echo 'failed to commit updates.';
    fi
  fi
  echo "--------------------"
done
