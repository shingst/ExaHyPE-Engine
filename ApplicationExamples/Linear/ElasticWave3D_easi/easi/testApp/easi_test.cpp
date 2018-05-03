#include <iostream>
#include "easi/YAMLParser.h"
#include "easi/ResultAdapter.h"
#include "asagi_reader.h"

struct ElasticMaterial{
  double rho, cs, cp ;
};

int main (){
  //getSingeCoordinate
  easi::Query query (2,3);
  query.x(0,0) = 1.0;
  query.x(0,1) = 2.0;
  query.x(0,2) = 3.0;
  query.x(1,0) = 3.0;
  query.x(1,1) = 4.0;
  query.x(1,2) = 5.0;


  AsagiReader asagiReader("");
  easi::YAMLParser parser(3,&asagiReader);
  //  std::string fileName=;
  easi::Component* model = parser.parse("101_asagi.yaml");

  ElasticMaterial material[2];
  easi::ArrayOfStructsAdapter<ElasticMaterial> adapter(material);
  adapter.addBindingPoint("lambda", &ElasticMaterial::rho);
  adapter.addBindingPoint("mu",     &ElasticMaterial::cs);
  adapter.addBindingPoint("rho",    &ElasticMaterial::cp);
  
  model->evaluate(query, adapter);
  
  delete model;
  
  for (unsigned j = 0; j < 2; ++j) {
    std::cout << material[j].cp << " " << material[j].rho << " " << material[j].cs << std::endl;
  }  
  
}
