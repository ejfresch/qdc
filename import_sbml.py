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


import sys
import re
from os.path import isfile
sys.path.append('./installer/libsbml/lib/python/dist-packages/libsbml')

from libsbml import *



def import_sbml(f_inp):

    #definition avogardo constant
    avogadro=float(602214179000000000000000)

    #I create a SBML Reader object
    sbml_reader=SBMLReader()
    
    #The function returns a SBMLDocument object
    sbml_doc=sbml_reader.readSBML(f_inp) 
    
    n_of_errors=sbml_doc.getNumErrors()
    if n_of_errors > 0:
        error=":: ERROR! This is not a SBML document!"
        return error
        
    version=sbml_doc.getVersion()
    level=sbml_doc.getLevel()
    if ((level !=2) or (version != 4)):
        print ":: WARNING: this software is designed for SBML level 2 version 4!"
    

    #The function returns a Model object
    model=sbml_doc.getModel()
    
     
    
    n_of_compartments=model.getNumCompartments()
    if n_of_compartments != 1:
        error=":: ERROR! Too many compartments!"
        return error

    for compartment in model.getListOfCompartments():
        volume=float(compartment.getVolume())
        my_volume="volume, "+str(volume)
    
    #Species and initial amounts
    list_of_species=[]
    initial_no_molec=[]
    null_species=0
    
    for species in model.getListOfSpecies():
        list_of_species.append(species.getId())
	if species.isSetInitialConcentration():
		#print float(volume)
		#print float(species.getInitialConcentration())
		#print avogadro
		#print float(species.getInitialConcentration())*float(volume)*avogadro
        	initialization_species="%s, 0.0, %d" % (species.getId(),int(float(species.getInitialConcentration())*float(volume)*avogadro))
		#print initialization_species
        	initial_no_molec.append(initialization_species)
    	elif species.isSetInitialAmount():
		initialization_species="%s, 0.0, %d" % (species.getId(),int(species.getInitialAmount()))
        	initial_no_molec.append(initialization_species)
    
    my_initialization="\n".join(initial_no_molec)


   
    #Global parameters
    dict_parameters={}

    for global_parameter in model.getListOfParameters():
        dict_parameters[global_parameter.getId()]=global_parameter.getValue()
    

    #Events
    dict_events={}
    list_of_events=[]
    list_of_events_init=[]
    dict_associations={}
    
    if model.getNumEvents() != 0:
        for idx_event, event in enumerate(model.getListOfEvents()):
            trigger=event.getTrigger()
            trig=formulaToString(trigger.getMath())
            pat=re.compile('\d+\.*\d*')
            results_0=re.findall(pat, trig)
            last_position=len(results_0)-1
            for assignment in range(event.getNumEventAssignments()):
                e_assignment=event.getEventAssignment(assignment)
                
                if dict_events.has_key(e_assignment.getVariable()):
                    var=str(e_assignment.getVariable())
                    list_of_events.append("1, "+str(float(results_0[last_position]))+", $"+var+", "+str(dict_parameters[formulaToString(e_assignment.getMath())]))
                    dict_associations[float(results_0[last_position])]=idx_event
                elif list_of_species.count(str(e_assignment.getVariable()))==1:
                    formula_assignment=formulaToString(e_assignment.getMath())
                    list_of_events.append("assegnaz")
                    pattern_component=re.compile('\+ (\w+)\)')
                    component=re.findall(pattern_component, formula_assignment)
                    my_initialization=my_initialization+"\n"+str(e_assignment.getVariable())+", "+str(float(results_0[last_position]))+", "+str(int(dict_parameters[component[0]]*avogadro))+"\n"


                else:
                    
                    var=str(e_assignment.getVariable())
                    list_of_events_init.append("0, $"+var+", "+str(dict_parameters[var]))
                    list_of_events.append("1, "+str(float(results_0[last_position]))+", $"+var+", "+str(dict_parameters[formulaToString(e_assignment.getMath())]))
                    dict_events[var]="True"      
                    dict_associations[float(results_0[last_position])]=idx_event
                    


        dict_list=dict_associations.keys()
        dict_list.sort()
        
        
        my_events="\n".join(list_of_events_init)
        
        
        for item_0 in dict_list:
            my_events=my_events+"\n"+list_of_events[dict_associations[item_0]]

    else:
        my_events=""


    #Reactions
    list_of_reactions=[]
    
    for idx_reaction, reaction in enumerate(model.getListOfReactions()):

        #Parameter
        if reaction.isSetFast() and reaction.getFast():
            parameter_value="-"
                            

        else:

            kinetic_law=reaction.getKineticLaw()
            if kinetic_law.getNumParameters() > 1:
                error=":: ERROR! Wrong number of parameters! (reaction %d)" % idx_reaction
                return error
            elif kinetic_law.getNumParameters() == 0:
                current_formula=kinetic_law.getFormula()
                times=re.compile(' \* ')
                results_1=re.split(times,current_formula)
                for item_1 in results_1:
                    if dict_parameters.has_key(item_1):
                        if dict_events.has_key(item_1):
                            parameter_value="$"+item_1
                        else:
                            parameter_value=str(dict_parameters[item_1])
            else:
                current_parameter=kinetic_law.getParameter(0)
                parameter_value=str(current_parameter.getValue())
                
        
        #Reactants
        
        list_of_reactants=[]

        
        if reaction.getNumReactants()==0:
            if null_species==0:
                null_species=1
                #list_of_species.append("null")
            list_of_reactants.append("null")
        
        else:
            for reactant in range(reaction.getNumReactants()):
        
		current_reactant=reaction.getReactant(reactant)
		stoichiometry=current_reactant.getStoichiometry()
                
                if reaction.isSetFast() and reaction.getFast():
                        item="%d %s" % (stoichiometry, current_reactant.getSpecies())    
		elif (float(reaction.getNumReactants())==1) and float(stoichiometry)== 2.0:
                    item="%d %s" % (stoichiometry, current_reactant.getSpecies())
                elif float(stoichiometry) != 1.0: 
			#print "qui"
			error= ":: ERROR! (stoichiometry)"
			return error
		else:
                    item=current_reactant.getSpecies()
                
                list_of_reactants.append(item)   
        


        #Products
        
        list_of_products=[]

        if reaction.getNumProducts()==0:
            if null_species==0:
                null_species=1
                #list_of_species.append("null")

            list_of_products.append("null")

        else:
            for product in range(reaction.getNumProducts()):
        	current_product=reaction.getProduct(product)
		stoichiometry=current_product.getStoichiometry()
		#if reaction.isSetFast() and reaction.getFast():
		if int(stoichiometry)>1:                
			item="%d %s" % (stoichiometry, current_product.getSpecies())    
		#elif float(stoichiometry) != 1.0: 			
		#	print "QUI"		        
		#	error= ":: ERROR! (stoichiometry)"
		#	return error
		else:
                    item=current_product.getSpecies()
		         
		
		list_of_products.append(item)   



        my_reaction=parameter_value+", "+" + ".join(list_of_reactants)+" > "+" + ".join(list_of_products)    
        list_of_reactions.append(my_reaction)
    

    my_reactions="\n".join(list_of_reactions)

    list_of_variables=[]
    for key in dict_events.keys():
        key="$"+key
        list_of_variables.append(key)
     
    my_variables=", ".join(list_of_variables)   
        
    my_species=", ".join(list_of_species)


    

    model_ok=my_species+"\n\n"+my_volume+"\n\n"+"time>1"+"\n\n"+my_reactions+"\n\n"+my_initialization+"\n\n"+my_variables+"\n\n"+my_events+"\n\n"
    return model_ok
    


if __name__ =='__main__':

    if len(sys.argv) != 2:
        print "usage: %s <file_input>" % sys.argv[0]
        sys.exit()


    file=sys.argv[1]
    if isfile(file):
        translation=import_sbml(file)
        print translation
    else:
        print ":: ERROR! This file does not exist!"


