#!/usr/bin/python                          

# Copyright 2009-2012 - Luca Freschi <l.freschi@gmail.com>                                                                       

# This file is part of QDC.                                                                                                
 
# QDC is free software: you can redistribute it and/or modify                                                               
# it under the terms of the GNU General Public License as published by                                                      
# the Free Software Foundation, either version 3 of the License, or                                                         
# (at your option) any later version.                                                                                       

# This program is distributed in the hope that it will be useful,                                                           
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                                            
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                             
# GNU General Public License for more details.                                                                              

# You should have received a copy of the GNU General Public License                                                         
# along with this program.  If not, see <http://www.gnu.org/licenses/>.   

import re
import sys
sys.path.append('./installer/libsbml/lib/python/dist-packages/libsbml')

from libsbml import *



def export_sbml(inp):
    
    #definition avogardo constant
    avogadro=float(602214179000000000000000)


    list_inp=[]
    list_inp=inp
    
    to_do=7

    counter_reaction=0
    counter_parameter=0
    counter_addition=0
    counter_event=0
    
    dict_species={}
    dict_parameters={}
    
    dict_associations={}


    section=to_do


    sbml_document=SBMLDocument()

    model=sbml_document.createModel("my_model")
    compartment=model.createCompartment()
    compartment.setId("comp")

    per_second=model.createUnitDefinition()
    per_second.setId("per_second")
    sec=per_second.createUnit()
    sec.setKind(UNIT_KIND_SECOND)
    sec.setExponent(-1)

    litre_per_mole_per_second=model.createUnitDefinition()
    litre_per_mole_per_second.setId("litre_per_mole_per_second")
    mol=litre_per_mole_per_second.createUnit()
    mol.setKind(UNIT_KIND_MOLE)
    mol.setExponent(-1)
    lit=litre_per_mole_per_second.createUnit()
    lit.setKind(UNIT_KIND_LITRE)
    lit.setExponent(1)
    sec_2=litre_per_mole_per_second.createUnit()
    sec_2.setKind(UNIT_KIND_SECOND)
    sec_2.setExponent(-1)

    mole_per_litre_per_second=model.createUnitDefinition()
    mole_per_litre_per_second.setId("mole_per_litre_per_second")
    mol_2=mole_per_litre_per_second.createUnit()
    mol_2.setKind(UNIT_KIND_MOLE)
    mol_2.setExponent(1)
    lit_2=mole_per_litre_per_second.createUnit()
    lit_2.setKind(UNIT_KIND_LITRE)
    lit_2.setExponent(-1)
    sec_3=mole_per_litre_per_second.createUnit()
    sec_3.setKind(UNIT_KIND_SECOND)
    sec_3.setExponent(-1)





    for idx_item_0, item_0 in enumerate(list_inp):
        
        
        item_0=item_0.replace(' ', '')
        item_0=item_0.lower()
        #print item_0
        #print section
        hash_0=re.compile('#')
        match_hash=re.findall(hash_0, item_0)
        if match_hash:
            list_hash=re.split(hash_0, item_0)
            if list_hash[0] == '':
                continue
            else:
                item_0=list_hash[0]

            
        if item_0 == '':
            
            if (to_do > 1): 
                to_do-=1
                section=to_do
                continue
            else:
                break


        if section == 7:
            sep=re.compile(',')
            list_of_species=re.split(sep, item_0)
            for item_1 in list_of_species:
                if item_1=="null":
                    null_species=1
		elif item_1=="":
		    skip=1
                else:
                    species=model.createSpecies()
                    species.setName(item_1)
                    species.setId(item_1)
                    dict_species[item_1]=species
                    species.setCompartment("comp")
                    species.setInitialConcentration(0)

        if section == 6:
            list_vol=re.split(sep, item_0)
            compartment.setVolume(float(list_vol[1]))



        if section == 4:


        
            re_first=re.compile('(.+),(\d*)(\w+)>(.+)')
            match_first=re.findall(re_first, item_0)

            re_second=re.compile('(.+),(\d*)(\w+)\+(\d*)(\w+)>(.+)')
            match_second=re.findall(re_second, item_0)


          
            if match_first:
                
                reaction=model.createReaction()
                counter_reaction+=1
                id_reaction="reaction"+str(counter_reaction)
                reaction.setId(id_reaction)
                reaction.setReversible(False)
            
    
                #Reagents

                if str(match_first[0][2])=="null":
                    null_reaction=1

                else:
                    reagent=reaction.createReactant()
                    reagent.setSpecies(match_first[0][2])
                        
                    if (str(match_first[0][1])!="" and str(match_first[0][0])=='-'):
                        reagent.setStoichiometry(float(match_first[0][1]))
                    
                    elif (str(match_first[0][1])=='2'):
                        reagent=reaction.createReactant()
                        reagent.setSpecies(match_first[0][2])
                #Products 
                pat_prod=re.compile('(\d*)(\w+)')
                match_prod=re.findall(pat_prod,match_first[0][3])
                while(len(match_prod)!=0):
                    sto,prod=match_prod.pop(0)
                    if prod=="null":
                        prod_null=1
                    else:
                        product=reaction.createProduct()
                        product.setSpecies(prod)
                        if sto != '':
                            product.setStoichiometry(float(sto))
                    

                #Local parameters and kinetic law
                if str(match_first[0][0]) == '-':
                    reaction.setFast("True")
                else:
                    
                    k=re.compile('\$(\w+)')
                    match_k=re.findall(k, match_first[0][0])
                    if match_k:
                        parameter=model.createParameter()
                        id_parameter=match_k[0]
                        dict_parameters[id_parameter]=parameter
                        parameter.setId(id_parameter)
                        kl=reaction.createKineticLaw()
                        if(str(match_first[0][2])=="null"):
                            parameter.setUnits("mole_per_litre_per_second")
                            formula_first=id_parameter+'*'+compartment.getId()
                        elif (str(match_first[0][1])=='2'):
                            parameter.setUnits("litre_per_mole_per_second")
                            formula_first=id_parameter+'*'+compartment.getId()+'*'+match_first[0][2]+'*'+match_first[0][2]
                        else:
                            formula_first=match_k[0]+'*'+compartment.getId()+'*'+match_first[0][2]
                            parameter.setUnits("per_second")
                        ast_node_first=parseFormula(formula_first)
                        kl.setMath(ast_node_first)
                        continue
                
                    parameter_2=model.createParameter()
                    parameter_2.setValue(float(match_first[0][0]))
                    counter_parameter+=1
                    id_parameter="parameter"+str(counter_parameter)
                    parameter_2.setId(id_parameter)
                    parameter_2.setConstant(True)
                    if str(match_first[0][2])=="null":
                        parameter_2.setUnits("mole_per_litre_per_second")
                    else:
                        parameter_2.setUnits("per_second")
                    kl=reaction.createKineticLaw()
                    if str(match_first[0][2])=="null":
                        formula_first=id_parameter+'*'+compartment.getId()
                    elif (str(match_first[0][1])=='2'):
                        parameter_2.setUnits("litre_per_mole_per_second")
                        formula_first=id_parameter+'*'+compartment.getId()+'*'+match_first[0][2]+'*'+match_first[0][2]
                    else:
                        formula_first=id_parameter+'*'+compartment.getId()+'*'+match_first[0][2] 
                    ast_node_first=parseFormula(formula_first)
                    kl.setMath(ast_node_first)
        
            
                        

            elif match_second:

                
                
                reaction=model.createReaction()
                counter_reaction+=1
                id_reaction="reaction"+str(counter_reaction)
                reaction.setId(id_reaction)
                reaction.setReversible(False)
            


                #Reagents

                if str(match_second[0][2])=="null" or str(match_second[0][4])=="null":
                    msg=":: ERROR! Invalid reaction type (please check: null)"
                    return msg
                reagent_1=reaction.createReactant()
                reagent_1.setSpecies(match_second[0][2])
                if match_second[0][1] != '':
                    reagent_1.setStoichiometry(float(match_second[0][1]))
                reagent_2=reaction.createReactant()
                reagent_2.setSpecies(match_second[0][4])
                if match_second[0][3] != '':
                    reagent_2.setStoichiometry(float(match_second[0][3]))
               
                #Products
                pat_prod=re.compile('(\d*)(\w+)')
                match_prod=re.findall(pat_prod,match_second[0][5])
                while(len(match_prod)!=0):
                    sto,prod=match_prod.pop(0)
                    if prod=="null":
                        prod_null=1
                    else:
                        product=reaction.createProduct()
                        product.setSpecies(prod)
                        if sto != '':
                            product.setStoichiometry(float(sto))
                                
                
                        
                #Local parameters and kinetic law
                if str(match_second[0][0]) == '-':
                    reaction.setFast("True")
                else:
                    k=re.compile('\$(\w+)')
                    match_k=re.findall(k, match_second[0][0])
                
                
                    if match_k:

                        parameter=model.createParameter()
                        id_parameter=match_k[0]
                        dict_parameters[id_parameter]=parameter
                        parameter.setId(id_parameter)
                        parameter.setUnits("litre_per_mole_per_second")
                        kl=reaction.createKineticLaw()
                                              
                        formula_second=match_k[0]+'*'+compartment.getId()+'*'+match_second[0][2]+'*'+match_second[0][4]
                        ast_node_second=parseFormula(formula_second)
                        kl.setMath(ast_node_second)
                        continue
                
                    parameter_2=model.createParameter()
                    parameter_2.setValue(float(match_second[0][0]))
                    counter_parameter+=1
                    id_parameter="parameter"+str(counter_parameter)
                    parameter_2.setId(id_parameter)
                    parameter_2.setConstant(True)
                    parameter_2.setUnits("litre_per_mole_per_second")
                    kl=reaction.createKineticLaw()
                    formula_second=id_parameter+'*'+compartment.getId()+'*'+match_second[0][2]+'*'+match_second[0][4] 
                    ast_node_second=parseFormula(formula_second)
                    kl.setMath(ast_node_second)

            
            else:
                
                print ":: ERROR! (reactions)"
                sys.exit()


        if section == 3:
            sep_init=re.compile(',')
            data_init=re.split(sep_init, item_0)
            if len(data_init)!=3:
                print ":: ERROR (Species init)"
                return
            if float(data_init[1])==float(0.0):
		#print data_init[2]
		#print float(data_init[2])/avogadro
                dict_species[data_init[0]].setInitialConcentration((float(data_init[2])/avogadro)/float(list_vol[1]))
            else:
                parameter=model.createParameter()
                counter_addition+=1
                id_parameter="a"+str(counter_addition)
                parameter.setId(id_parameter)
                parameter.setValue(float(data_init[2])/avogadro)
                parameter.setConstant(True)
                parameter.setUnits("mole")
                
                #Additions
                event=model.createEvent()
                counter_event+=1
                id_event="event"+str(counter_event)
                event.setId(id_event)
                formula="eq("+"t"+","+data_init[1]+")"
                ast_node_trigger=parseFormula(formula)
                node=ast_node_trigger.getChild(0)
                node.setType(AST_NAME_TIME)
                my_trigger=Trigger(2, 4)
                my_trigger.setMath(ast_node_trigger)
                event.setTrigger(my_trigger)
                assignment=event.createEventAssignment()
                formula_assignment='('+data_init[0]+'* comp +'+id_parameter+') / comp'
                ast_node_assignment=parseFormula(formula_assignment)
                assignment.setMath(ast_node_assignment)
                assignment.setVariable(data_init[0])
                

        if section == 1:
           	
            re_init_var=re.compile('^0,\$(\w+),(\d+\.?\d*)$')
            re_ev_var=re.compile('^1,(\d+\.?\d*),\$(\w+),(\d+\.?\d*)$')

            match_init_var=re.findall(re_init_var, item_0)
            match_ev_var=re.findall(re_ev_var, item_0)

            if (len(match_init_var)==1 and len(match_init_var[0])==2):
                if dict_parameters.has_key(match_init_var[0][0]):
                    parameter=dict_parameters[match_init_var[0][0]]
                    parameter.setValue(float(match_init_var[0][1]))
                    parameter.setConstant(False)
            
            
            elif (len(match_ev_var) == 1 and len(match_ev_var[0])==3):               
		parameter=model.createParameter()
                counter_parameter+=1
                id_parameter="parameter"+str(counter_parameter)
                dict_parameters[id_parameter]=parameter
                parameter.setId(id_parameter)
                parameter.setValue(float(match_ev_var[0][2]))
                event=model.createEvent()
                counter_event+=1
                id_event="event"+str(counter_event)
                event.setId(id_event)
                formula="geq("+"t"+","+match_ev_var[0][0]+")"
                ast_node_trigger=parseFormula(formula)
                node=ast_node_trigger.getChild(0)
                node.setType(AST_NAME_TIME)
                my_trigger=Trigger(2, 4)
                my_trigger.setMath(ast_node_trigger)
                event.setTrigger(my_trigger)
                assignment=event.createEventAssignment()
                formula_assignment=id_parameter
                ast_node_assignment=parseFormula(formula_assignment)                
		assignment.setMath(ast_node_assignment)
                assignment.setVariable(match_ev_var[0][1])
                dict_associations[id_parameter]=match_ev_var[0][1]
                
            else:
                print ":: ERROR! (events)"
                return
        


    for key in dict_associations.keys():
        id_u=dict_associations[key]
        parameter_u=dict_parameters[id_u]
        unit=parameter_u.getUnits()
        parameter=dict_parameters[key]
        parameter.setUnits(unit)


    model_ok=writeSBMLToString(sbml_document)
    return model_ok

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "usage: %s <file_input>" % sys.argv[0]
        sys.exit()

    file=sys.argv[1]
    handler=open(file, 'r')
    my_list=[]
    for line in handler:
        line=line.replace('\n','')
        my_list.append(line)
        
    handler.close()
    stri=export_sbml(my_list)
    print stri
