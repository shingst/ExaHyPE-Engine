package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;

public class UserDefinedADER_DGinFortran implements Solver {
  private String _projectName;
  private String _solverName;
 
  public UserDefinedADER_DGinFortran(String projectName, String solverName) {
    _projectName = projectName;
    _solverName = solverName;
  }
  
    
  @Override
  public String getSolverName() {
    return _solverName;
  }
  
  @Override
  public void writeAbstractHeader(BufferedWriter writer) throws IOException {
    // TODO Auto-generated method stub
    
  }
  
  @Override
  public void writeAbstractImplementation(BufferedWriter writer) throws IOException {
    // TODO Auto-generated method stub
    
  }
  
  public void writeHeader(java.io.BufferedWriter writer)
      throws java.io.IOException {
    // @todo
    System.err.println("not implemented yet\n");
  }

  public void writeUserImplementation(java.io.BufferedWriter writer) throws java.io.IOException {
    // @todo Implement
    System.err.println("not implemented yet\n");
  }
  public void writeUserPDE(java.io.BufferedWriter writer)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("not implemented yet\n");
  }
  public void writeTypesDef(java.io.BufferedWriter writer)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("not implemented yet\n");
  }
  
  @Override
  public boolean supportsVariables() {
    return false;
  }
}
