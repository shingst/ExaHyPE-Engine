#!/usr/bin/env groovy
util = load "Miscellaneous/Jenkins/Util.groovy"

def subject(){
    commit_nr = GIT_COMMIT.take(7)
    return '''ExaHyPE-Engine | Pipeline has failed for '''+BRANCH_NAME+''' | '''+commit_nr
}

def body(){
    url=GIT_URL.replaceAll('git@','https://').replaceAll('gitlab.lrz.de:','gitlab.lrz.de/').replaceAll('ExaHyPE-Engine.git','ExaHyPE-Engine')
    commit_nr=GIT_COMMIT
    commit_url=url+"/commit/"+commit_nr
    commit_nr_small=GIT_COMMIT.take(7)
    
    def body='''Your pipeline has failed.\n Project: ExaHyPE-Engine '''+url+'''\n Branch: '''+BRANCH_NAME+''' \n Commit: '''+commit_nr_small+'''  '''+commit_url+'''\n See pipeline: '''+RUN_DISPLAY_URL

    return body
}

return this
