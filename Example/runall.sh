# /bin/bash
pwd1=`pwd`
cd $pwd1 ; cd ..
pwd2=`pwd`
cd $pwd1/Plot ; $pwd2/Bin/amesp CH4.aip ; $pwd2/Bin/amesp NaS.aip ; $pwd2/Bin/amesp NH3.aip ; $pwd2/Bin/amesp H12.aip ; rm *.mo
cd $pwd1/HF/RHF ; $pwd2/Bin/amesp C6H6.aip ; $pwd2/Bin/amesp H60.aip
cd $pwd1/HF/ROHF ; $pwd2/Bin/amesp *.aip
cd $pwd1/HF/UHF ; $pwd2/Bin/amesp *.aip
cd $pwd1/CI/CIS ; $pwd2/Bin/amesp *.aip
cd $pwd1/CI/CID ; $pwd2/Bin/amesp *.aip
cd $pwd1/CI/CISD ; $pwd2/Bin/amesp *.aip
cd $pwd1/DFT/B3LYP ; $pwd2/Bin/amesp HF.aip ; $pwd2/Bin/amesp H2O.aip ; $pwd2/Bin/amesp HO.aip
cd $pwd1/DFT/B3PW91 ; $pwd2/Bin/amesp HF.aip ; $pwd2/Bin/amesp H2O.aip ; $pwd2/Bin/amesp HO.aip
cd $pwd1/DFT/PW91 ; $pwd2/Bin/amesp HF.aip ; $pwd2/Bin/amesp H2O.aip ; $pwd2/Bin/amesp HO.aip
cd $pwd1/DFT/BLYP ; $pwd2/Bin/amesp HF.aip ; $pwd2/Bin/amesp H2O.aip ; $pwd2/Bin/amesp HO.aip
cd $pwd1/DFT/PBE ; $pwd2/Bin/amesp HF.aip ; $pwd2/Bin/amesp H2O.aip ; $pwd2/Bin/amesp HO.aip
cd $pwd1/DFT/PBE0 ; $pwd2/Bin/amesp HF.aip ; $pwd2/Bin/amesp H2O.aip ; $pwd2/Bin/amesp HO.aip
cd $pwd1/DFT/XAlpha ; $pwd2/Bin/amesp HF.aip ; $pwd2/Bin/amesp H2O.aip ; $pwd2/Bin/amesp HO.aip
cd $pwd1/DFT/LSDA ; $pwd2/Bin/amesp HF.aip ; $pwd2/Bin/amesp H2O.aip ; $pwd2/Bin/amesp HO.aip
cd $pwd1/MBPT/MP2 ; $pwd2/Bin/amesp C2H6.aip ; $pwd2/Bin/amesp COFH.aip
cd $pwd1/MBPT/MP3 ; $pwd2/Bin/amesp C2H4.aip ; $pwd2/Bin/amesp H13.aip
cd $pwd1/TD/TDHF ; $pwd2/Bin/amesp CHOH.aip
