#!/usr/bin/env groovy
import jenkins.model.Jenkins

pipeline {
    agent { label 'Linux-Cluster' }
    triggers {
        pollSCM('H H * * *')
     	//cron('H H * * *')
    }
    stages {
	stage ('Checkout') {
	    steps {
		sh 'echo Running as $(id -un)@$(hostname -f)'
		checkout scm
		sh 'pwd'
		sh 'ls'

		script {
		config = load "Miscellaneous/Jenkins/Config.groovy"
		matrix = load "Miscellaneous/Jenkins/Matrix.groovy"
		run = load "Miscellaneous/Jenkins/Run.groovy"
		build = load "Miscellaneous/Jenkins/Build.groovy"
		}
		sh util.getModuleCode() + 'Submodules/updateSubmodules.sh -s'
		sh util.getModuleCode() + 'Submodules/updateSubmodules.sh'

	    }
	}

	stage ('Build') {
	    steps {
		sh 'pwd'
		sh 'ls'
		sh 'ls ApplicationExamples || true'
		// Gather all test cases.
		script {
		    testCases = matrix.getMatrix()
		    def buildMatrix = build.assembleBuildMatrix(testCases, "${workspace}")
		    parallel buildMatrix
		}
	    }
	}
	stage ('Run') {
	    steps {
		script {
		    def runMatrix = run.assembleRunMatrix(testCases, "${workspace}")
		    parallel runMatrix
		}
	    }
	}
        stage ('Deploy'){
          when { branch "jenkins" }     
            steps{
              sh "git checkout jenkins"
              sh "git pull"
              sh "git checkout release"
              sh "git pull"
              sh "git merge jenkins"
//              sh 'git commit -ma "Jenkins: Deployed master to release branch"'
              sh "git push origin release"
            }
          }
        }
    post {
	always {
	    cleanWs()
	}
	success {
	    updateGitlabCommitStatus(name: BRANCH_NAME, state: 'success')
	}
	failure {
	    updateGitlabCommitStatus(name: BRANCH_NAME, state: 'failed')
	}
    }
    options {
	gitLabConnection('Exahype-Gitlab')
    }
}
