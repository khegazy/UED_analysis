#include "/reg/neh/home/khegazy/baseScripts/atomClass.h"
#include "/reg/neh/home/khegazy/baseScripts/moleculeClass.h"
#include "/reg/neh/home/khegazy/baseScripts/molEnsembleMC.h"
#include "/reg/neh/home/khegazy/baseScripts/diffractionClass.h"
#include "/reg/neh/home/khegazy/baseScripts/plotClass.h"
#include "/reg/neh/home/khegazy/baseScripts/imageProcessing.h"
#include "/reg/neh/home/khegazy/baseScripts/saveClass.h"


enum radicalEnums {nitrobenzene, phenoxyRadical, phenylRadical};


//////////////////////////
/////  Nitrobenzene  /////
//////////////////////////

class NBZMCclass : public MOLENSEMBLEMCclass {

    public:
        NBZMCclass(long int seed, string pdfPath, string pdfNames) 
                : MOLENSEMBLEMCclass(seed, pdfPath, pdfNames) {}
        NBZMCclass(long int seed) : MOLENSEMBLEMCclass(seed) {}

        void buildMolecule(MOLECULEclass &molecule, 
            std::map<std::string, double> inpVals) {

	  TVector3 position;

	  position.SetXYZ(0.1102854501, -0.0000017476, -0.0960606251); 
          molecule.addAtom(new ATOMclass("C0", C, angs_to_au*position));

          position.SetXYZ(-0.0233994016, 0.0000115614, 1.2949218659);
          molecule.addAtom(new ATOMclass("C1", C, angs_to_au*position));

          position.SetXYZ(1.1339053844, -0.0000132121, 2.0761968001);
          molecule.addAtom(new ATOMclass("C2", C, angs_to_au*position));
          
          position.SetXYZ(2.3959032389, 0.0000068205, 1.4661586557);
          molecule.addAtom(new ATOMclass("C3", C, angs_to_au*position));
          
          position.SetXYZ(2.5060086613, 0.0000130021, 0.0687992624);
          molecule.addAtom(new ATOMclass("C4", C, angs_to_au*position));
          
          position.SetXYZ(1.3577715577, -0.0000194139, -0.7257476393);
          molecule.addAtom(new ATOMclass("C5", C, angs_to_au*position));
          
          position.SetXYZ(-1.0115675324, -0.0000031181, 1.7353472487);
          molecule.addAtom(new ATOMclass("H1", H, angs_to_au*position));
          
          position.SetXYZ(1.0509346293, 0.0000028455, 3.1572922713);
          molecule.addAtom(new ATOMclass("H2", H, angs_to_au*position));
          
          position.SetXYZ(3.2916307024, 0.0000046106, 2.0783087949);
          molecule.addAtom(new ATOMclass("H3", H, angs_to_au*position));
          
          position.SetXYZ(3.4832038999, -0.0000202132, -0.4010426025);
          molecule.addAtom(new ATOMclass("H4", H, angs_to_au*position));
          
          position.SetXYZ(1.4096531886, -0.0000197767, -1.8063691645);
          molecule.addAtom(new ATOMclass("H5", H, angs_to_au*position));
          
          position.SetXYZ(-1.1003992842, 0.0000026617, -0.9236463600);
          molecule.addAtom(new ATOMclass("N0", N, angs_to_au*position));
          
          position.SetXYZ(-2.2217452071, -0.0000033365, -0.3388917257);
          molecule.addAtom(new ATOMclass("O1", O, angs_to_au*position));
          
          position.SetXYZ(-0.9626768830, 0.0000054957, -2.1807621240);
          molecule.addAtom(new ATOMclass("O2", O, angs_to_au*position));
       }
};





/////////////////////////////
/////  Phenoxy Radical  /////
/////////////////////////////

class PNOXYRadMCclass : public MOLENSEMBLEMCclass {

    public:
        PNOXYRadMCclass(long int seed, string pdfPath, string pdfNames) 
                : MOLENSEMBLEMCclass(seed, pdfPath, pdfNames) {}
        PNOXYRadMCclass(long int seed) : MOLENSEMBLEMCclass(seed) {}

        void buildMolecule(MOLECULEclass &molecule, 
            std::map<std::string, double> inpVals) {

          TVector3 position;

          position.SetXYZ(0.0032109344, -0.0000004058, -0.0898762103);
          molecule.addAtom(new ATOMclass("C0", C, angs_to_au*position));

          position.SetXYZ(-0.0216886180, 0.0000018685, 1.3505315895);
          molecule.addAtom(new ATOMclass("C1", C, angs_to_au*position));

          position.SetXYZ(1.1581806546, -0.0000064309, 2.0759959342);
          molecule.addAtom(new ATOMclass("C2", C, angs_to_au*position));

          position.SetXYZ(2.3992790292, 0.0000045587, 1.4069028434);
          molecule.addAtom(new ATOMclass("C3", C, angs_to_au*position));

          position.SetXYZ(2.4562711928, 0.0000044575, -0.0019018439);
          molecule.addAtom(new ATOMclass("C4", C, angs_to_au*position));

          position.SetXYZ(1.2868045572, -0.0000083038, -0.7440087670);
          molecule.addAtom(new ATOMclass("C5", C, angs_to_au*position));

          position.SetXYZ(-0.9925490793, 0.0000029311, 1.8332329709);
          molecule.addAtom(new ATOMclass("H1", H, angs_to_au*position));

          position.SetXYZ(1.1338536441, -0.0000061065, 3.1608456117);
          molecule.addAtom(new ATOMclass("H2", H, angs_to_au*position));

          position.SetXYZ(3.3196837485, 0.0000076177, 1.9818618723);
          molecule.addAtom(new ATOMclass("H3", H, angs_to_au*position));

          position.SetXYZ(3.4205708327, 0.0000042510, -0.4995172989);
          molecule.addAtom(new ATOMclass("H4", H, angs_to_au*position));

          position.SetXYZ(1.2946410930, -0.0000117188, -1.8282061290);
          molecule.addAtom(new ATOMclass("H5", H, angs_to_au*position));

          position.SetXYZ(-1.0974837183, 0.0000033834, -0.7778726589);
          molecule.addAtom(new ATOMclass("O0", O, angs_to_au*position));
       }
};


/////  NO  /////

class NOMCclass : public MOLENSEMBLEMCclass {

    public:
        NOMCclass(long int seed, string pdfPath, string pdfNames) 
                : MOLENSEMBLEMCclass(seed, pdfPath, pdfNames) {}
        NOMCclass(long int seed) : MOLENSEMBLEMCclass(seed) {}

        void buildMolecule(MOLECULEclass &molecule, 
            std::map<std::string, double> inpVals) {

          TVector3 position;

          double rad = 1.15;

          position.SetXYZ(0, 0, 0);
          molecule.addAtom(new ATOMclass("N", N, angs_to_au*position));

          position.SetXYZ(0, 0, rad);
          molecule.addAtom(new ATOMclass("O", O, angs_to_au*position));
       }
};





////////////////////////////
/////  Phenyl Radical  /////
////////////////////////////

class PNLRadMCclass : public MOLENSEMBLEMCclass {

    public:
        PNLRadMCclass(long int seed, string pdfPath, string pdfNames) 
                : MOLENSEMBLEMCclass(seed, pdfPath, pdfNames) {}
        PNLRadMCclass(long int seed) : MOLENSEMBLEMCclass(seed) {}

        void buildMolecule(MOLECULEclass &molecule, 
            std::map<std::string, double> inpVals) {

          TVector3 position;

          position.SetXYZ(0.0454389558, 0.0000023054, 0.0262344660);
          molecule.addAtom(new ATOMclass("C0", C, angs_to_au*position));

          position.SetXYZ(-0.0201209563, 0.0000049732, 1.4080616360);
          molecule.addAtom(new ATOMclass("C1", C, angs_to_au*position));

          position.SetXYZ(1.2055307189, -0.0000071576, 2.1025828445);
          molecule.addAtom(new ATOMclass("C2", C, angs_to_au*position));

          position.SetXYZ(2.4152626012, 0.0000027254, 1.3944561144);
          molecule.addAtom(new ATOMclass("C3", C, angs_to_au*position));

          position.SetXYZ(2.4236459857, 0.0000038740, -0.0072672627);
          molecule.addAtom(new ATOMclass("C4", C, angs_to_au*position));

          position.SetXYZ(1.2093525254, -0.0000064298, -0.7214750904);
          molecule.addAtom(new ATOMclass("C5", C, angs_to_au*position));

          position.SetXYZ(-0.9644948488, 0.0000039457, 1.9422592979);
          molecule.addAtom(new ATOMclass("H1", H, angs_to_au*position));

          position.SetXYZ(1.2081465319, -0.0000089812, 3.1886885322);
          molecule.addAtom(new ATOMclass("H2", H, angs_to_au*position));

          position.SetXYZ(3.3552065702, 0.0000033593, 1.9371854151);
          molecule.addAtom(new ATOMclass("H3", H, angs_to_au*position));

          position.SetXYZ(3.3655343316, 0.0000049005, -0.5480769877);
          molecule.addAtom(new ATOMclass("H4", H, angs_to_au*position));

          position.SetXYZ(1.1997923270, -0.0000066834, -1.8064271854);
          molecule.addAtom(new ATOMclass("H5", H, angs_to_au*position));
       }
};


/////  NO2  /////

class NO2MCclass : public MOLENSEMBLEMCclass {

    public:
        NO2MCclass(long int seed, string pdfPath, string pdfNames) 
                : MOLENSEMBLEMCclass(seed, pdfPath, pdfNames) {}
        NO2MCclass(long int seed) : MOLENSEMBLEMCclass(seed) {}

        void buildMolecule(MOLECULEclass &molecule, 
            std::map<std::string, double> inpVals) {

          TVector3 position;

          double rad = 1.197;
          double ang = 134.3*PI/180;

          position.SetXYZ(0, 0, 0);
          molecule.addAtom(new ATOMclass("N0", N, angs_to_au*position));

          position.SetXYZ(rad*cos(ang/2), 0, rad*sin(ang/2));
          molecule.addAtom(new ATOMclass("O1", O, angs_to_au*position));

          position.SetXYZ(rad*cos(ang/2), 0, -rad*sin(ang/2));
          molecule.addAtom(new ATOMclass("O2", O, angs_to_au*position));
       }
};

