package eu.exahype.solvers;

import java.util.Set;
import eu.exahype.kernel.FiniteVolumesKernel;

/**
 * @author Dominic E. Charrier
 */
public class GenericFiniteVolumesGodunovInC extends GenericFiniteVolumesInC {
    public GenericFiniteVolumesGodunovInC(
            String projectName, String solverName, int dimensions, int numberOfVariables, int numberOfParameters, Set<String> namingSchemeNames, int patchSize, boolean enableProfiler, boolean hasConstants, FiniteVolumesKernel kernel) {
        super("godunov", projectName, solverName, dimensions,numberOfVariables,numberOfParameters,namingSchemeNames,patchSize,enableProfiler,hasConstants,kernel);
    }
}

