#!/usr/bin/env groovy
import jenkins.model.Jenkins

pipeline {
    agent { label 'Linux-Cluster' }
    // triggers {
    // 	cron('H H * * *')
    // }
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
		sh util.getModuleCode() + 'Peano/updatePeano.sh -s'
		sh util.getModuleCode() + 'Peano/updatePeano.sh'

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

    }

    post {
	always {
	    cleanWs()
	}
	success {
	    updateGitlabCommitStatus(name: 'jenkins', state: 'success')
	}
	failure {
	    updateGitlabCommitStatus(name: 'jenkins', state: 'failed')
	}
    }
    options {
	gitLabConnection('Exahype-Gitlab')
    }
}
