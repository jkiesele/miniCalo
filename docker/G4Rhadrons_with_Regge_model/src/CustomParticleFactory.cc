#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <iterator>
#include <map>
#include <set>

#include "CustomPDGParser.h"
#include "CustomParticle.h"
#include "CustomParticleFactory.h"
#include "G4DecayTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"

using namespace CLHEP;

bool CustomParticleFactory::loaded = false;
std::set<G4ParticleDefinition *> CustomParticleFactory::m_particles;

bool CustomParticleFactory::isCustomParticle(G4ParticleDefinition *particle)
{
  return (m_particles.find(particle)!=m_particles.end());
}

void CustomParticleFactory::loadCustomParticles()
{
  if(loaded) return;
  loaded = true;
  std::ifstream configFile("particles.txt");
  //  std::ifstream configFile("stophadrons.txt");
  G4String pType="custom";
  G4String pSubType="";
  G4double spectatormass;
  G4ParticleDefinition* spectator; 
  bool first = true;
  double mass;
  int pdgCode;
  std::string name,line;
  // This should be compatible IMO to SLHA 
  while(getline(configFile,line))
    {
      std::string::size_type beg_idx,end_idx;
      
      beg_idx = line.find_first_not_of("\t #");
      if(beg_idx > 0 && line[beg_idx-1] == '#') continue; 
      end_idx = line.find_first_of( "\t ", beg_idx);
      if (end_idx == std::string::npos) continue;
      pdgCode  = atoi(line.substr( beg_idx, end_idx - beg_idx ).c_str());
      
      beg_idx = line.find_first_not_of("\t ",end_idx);
      end_idx = line.find_first_of( "\t #", beg_idx);
      if (end_idx == std::string::npos) continue;
      mass  = atof(line.substr( beg_idx, end_idx - beg_idx ).c_str());
      
      beg_idx = line.find_first_not_of("\t# ",end_idx);
      end_idx = line.length();
      name = line.substr( beg_idx, end_idx - beg_idx );
      while(name.c_str()[0] == ' ') name.erase(0,1);
      while(name[name.size()-1] == ' ') name.erase(name.size()-1,1);
      
      //      std::cout << name << std::endl;
      
      if(abs(pdgCode) / 1000000 == 0)
	{
	  std::cout << "Pdg code too low " << pdgCode << " "<<abs(pdgCode) / 1000000  << std::endl;
	  continue;
	}
      
      pType="custom";
      if(CustomPDGParser::s_isRHadron(pdgCode)) pType = "rhadron";
      if(CustomPDGParser::s_isSLepton(pdgCode)) pType = "sLepton";
      if(CustomPDGParser::s_isMesonino(pdgCode)) pType = "mesonino";
      if(CustomPDGParser::s_isSbaryon(pdgCode)) pType = "sbaryon";

      std::cout<<"pType: "<<pType<<std::endl;
      std::cout<<"Charge of "<<name<<" is "<<CustomPDGParser::s_charge(pdgCode)<<std::endl;
      //      std::cout<<"s_isRMeson(pdgCode): "<<CustomPDGParser::s_isRMeson(pdgCode)<<std::endl;
      //      std::cout<<"s_isRBaryon(pdgCode): "<<CustomPDGParser::s_isRBaryon(pdgCode)<<std::endl;
      if(CustomPDGParser::s_isRMeson(pdgCode)) std::cout<<"Is R-meson"<<std::endl;
      if(CustomPDGParser::s_isRBaryon(pdgCode)) std::cout<<"Is R-baryon"<<std::endl;

      CustomParticle *particle  = new CustomParticle(
						     name,           mass * GeV ,        0.0*MeV,       eplus* CustomPDGParser::s_charge(pdgCode), 
						     (int)CustomPDGParser::s_spin(pdgCode)-1,              +1,             0,          
						     0,              0,             0,             
						     pType,               0,            +1, pdgCode,
						     true,            -1.0,          NULL );
      if (pType=="custom") {
	spectatormass = mass;
	spectator=particle;
	particle->SetCloud(0);
	particle->SetSpectator(0);
      }
      if(first){
	first = false;
	spectatormass = mass;
	spectator=particle;
	particle->SetCloud(0);
	particle->SetSpectator(0);
      } else {
	G4String cloudname = name+"cloud";
	G4String cloudtype = pType+"cloud";
	G4double cloudmass = mass-spectatormass;
	CustomParticle *tmpParticle  = new CustomParticle(
							  cloudname,           cloudmass * GeV ,        0.0*MeV,  0 , 
							  0,              +1,             0,          
							  0,              0,             0,             
							  cloudtype,               0,            +1, 0,
							  true,            -1.0,          NULL );
	particle->SetCloud(tmpParticle);
	particle->SetSpectator(spectator);

	G4cout<<name<<" being assigned "<<particle->GetCloud()->GetParticleName()<<" and "<<particle->GetSpectator()->GetParticleName()<<G4endl;
	G4cout<<"Masses: "
	      <<particle->GetPDGMass()/GeV<<" Gev, "
	      <<particle->GetCloud()->GetPDGMass()/GeV<<" GeV and "
	      <<particle->GetSpectator()->GetPDGMass()/GeV<<" GeV."
	      <<G4endl;
      }


      m_particles.insert(particle);
    }
  configFile.close();

  // Reading decays from file
  std::vector<std::vector<std::string>* > decays;
  std::ifstream decayFile("decays.txt");
  while(getline(decayFile,line))
    {
      std::istringstream is(line);
      std::vector<std::string>* txtvec = new std::vector<std::string>(std::istream_iterator<std::string>(is),std::istream_iterator<std::string>());
      decays.push_back(txtvec);
    }
  decayFile.close();

  // Looping over custom particles to add decays
  for (std::set<G4ParticleDefinition *>::iterator part=m_particles.begin();part!=m_particles.end();part++)
    {
      name=(*part)->GetParticleName();
      std::vector<std::vector<std::string> > mydecays;
      for (int i = 0; i!= decays.size(); i++)
	if(name==(*(decays[i]))[0]) // Is this decay for me?
	  mydecays.push_back(*(decays[i]));

      if (mydecays.size()>0) // Did I get any decays?
	{
	  int ndec=mydecays.size();
	  G4DecayTable* table = new G4DecayTable();
	  G4VDecayChannel** mode = new G4VDecayChannel*[ndec];
	  for (int i=0;i!=ndec;i++)
	    {
	      std::vector<std::string> thisdec=mydecays[i];
	      double branch = atof(thisdec[thisdec.size()-1].c_str()); // Reading branching ratio
	      thisdec.pop_back(); // Removing the number from the vector
	      if(thisdec.size()==3)
		mode[i] = new G4PhaseSpaceDecayChannel(thisdec[0],branch,2,thisdec[1],thisdec[2]);
	      else if(thisdec.size()==4)
		mode[i] = new G4PhaseSpaceDecayChannel(thisdec[0],branch,3,thisdec[1],thisdec[2],thisdec[3]);
	      else
		std::cout<<"Decay chain invalid!!!"<<std::endl;
	    }
	  for (G4int index=0; index <ndec; index++ ) table->Insert(mode[index]);
	  delete [] mode;
	  //	  (*part)->SetPDGLifeTime(0.1*ns);
	  //	  (*part)->SetPDGStable(false);
    	  (*part)->SetDecayTable(table);

	}
    }
  return;    
}

