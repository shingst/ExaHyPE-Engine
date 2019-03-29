#!/usr/bin/env groovy
import hudson.FilePath
import Config

def getMatrix() {
    // Warning: This might break down if more than one exahpye file is in a directory!
    def files = findFiles(glob: 'ApplicationExamples/TestCases/*/*.exahype*')
    files.each { d -> echo d.getPath() }
    testCases = []
    for (file in files) {
	testCases << new Config(exahypeFile: file.path,
				name: file.name)
    }
    return testCases
}

return this
