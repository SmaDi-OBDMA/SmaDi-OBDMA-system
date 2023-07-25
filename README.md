# Ontop-based OBDMA-sytem

## Introduction

This is the implementation of an OBDMA-system based on the OBDA-system-implementation ontop and a PostgresQL-database.
It shows the opportunities arising from modeling information in an OBDMA-system and is not meant to be a complete solution.

The system contains the ontology, a SPARQL-endpoint for stating individual queries and an user-interface for standard queries.

## Installation

The system is based on docker, therefore a docker installation is needed (https://www.docker.com/products/docker-desktop/) .
After starting docker, the programm can be executed by clicking on run.bat (Windows)  or executing run.sh (Mac OS) (via Terminal by executing "./run.sh").

(Starting up the system could take several minutes)


### Add an user to the docker-group

To run docker without administrator rights, the user have to be added to the docker-group. The following instructions are needed (for windows as administrator):

- search vor "computer management"

- System -> local users and groups -> groups -> docker-users

- click "add", then type the name of the account in the field "enter the object names to be used". After, click "check name" and "OK"

- restart your computer.


## SPARQL-templates corresponding to the presented use cases


### Maximum blocking stress

Natural query: Give me the maximum blocking stress with the information about the material class, material name, parameter name, value, unit, derivation, input parameter (if calculated) and a description if avaible.

#### SPARQL-query:

```
PREFIX : <urn:absolute/smadiont#>

SELECT ?matClass ?specimen ?parameter ?value ?unit ?derivation  ?input ?input_value ?input_unit ?description
WHERE { ?param  a  :maximum_blocking_stress ; 	
                :has_name ?parameter ;
                :has_value ?value ;
                :has_unit ?unit .
                         
        ?p      :has_parameter ?param ;
                :has_name ?specimen .   
        values  (?type ?matClass) {
                (:MSMA "MSMA")
                (:DE "DE")
                (:SMA "SMA")
                (:PC "PC")}
        ?p      a ?type.

optional{?param :determined_with ?m.		
        ?m      :has_name ?derivation.
optional{?m     :has_input ?e .	
        ?e      :has_name ?input.
optional{?e     :has_value ?input_value ;
                :has_unit ?input_unit .}}
optional{?m :has_description ?description .}}
}
```

### SMA - stress-strain curve depending on martensite content

Natural query: Give me the stress-strain curve depending on the martensite content with the information about the parameter names, unit, values and the martensite content as the condition. Choose the strain as the first parameter.

#### SPARQL-query:

```
PREFIX : <urn:absolute/smadiont#>

SELECT  ?specimen ?parameter1 ?value1 ?unit1 ?parameter2  ?value2 ?unit2   ?condition ?condition_unit ?condition_value
WHERE { ?a      a :stress-strain_curve_martensite.	 
  		?a      :has_tuple ?x .		
		?x      :curve_param1 ?param1.	
		?x      :curve_param2 ?param2.
		?param1 :has_name ?parameter1;
		        :has_value ?value1;
		        :has_unit ?unit1 .
		?param2 :has_name ?parameter2;
		        :has_value ?value2 ;
		        :has_unit ?unit2 .
		
  		?mat    :has_parameter ?param1 .
		?mat    :has_name ?specimen .
		
        ?a      :has_boundary_condition ?cond .
        ?cond   :has_name ?condition;
                :has_unit ?condition_unit;
                :has_value ?condition_value.
  
  	FILTER(?parameter1 = 'strain').	
}
```


### MSMA - maximum magnetic strain

Natural query: Give me the maximum magnetic strain with the information about the parameter name, unit, values and the derivation.


#### SPARQL-query:

```
PREFIX : <urn:absolute/smadiont#>

SELECT ?specimen ?parameter ?value ?unit ?derivation ?equ ?input ?input_value ?input_unit
WHERE { ?param  a  :maximum_strain_mag ; 	
                :has_name ?parameter ;
                :has_value ?value ;
                :has_unit ?unit .
                           
        ?p      :has_parameter ?param ;
                :has_name ?specimen .
           
optional{?param :determined_with ?m.		
        ?m      :has_name ?derivation;
				:has_equation ?equ.
optional{?m     :has_input ?e.
        ?e      :has_name ?input.
      	?e      :has_value ?input_value ;
      	        :has_unit ?input_unit .}}
} 
```


### PC - piezoelectric stress constant

Natural query: Give me g31 with information about the parameter name, unit, value and the derivation.

#### SPARQL-query:

```
PREFIX : <urn:absolute/smadiont#>

SELECT ?specimen ?paraname ?value ?unit ?derivation ?h ?descrip  

WHERE { ?param  :has_symbol "g";
                :has_row "3";
                :has_column "1";
	            :has_name ?paraname;
	            :has_unit ?unit; 
                :has_value ?value.            

        ?p      :has_parameter ?param ;
                :has_name ?specimen .  
optional{?param :determined_with ?h.
        ?h      :has_name ?derivation.}
}
```

### DE - elastic parameters

Natural query: Give me the elastic parameters for DE-Materials with the information about the parameter name, unit, value and the derivation.


#### SPARQL-query:

```
PREFIX : <urn:absolute/smadiont#>

SELECT ?specimen ?parameter ?value ?unit ?derivation  ?input ?input_value ?input_unit 
WHERE { ?param  a  :elastic_parameters  ; 	
                :has_name ?parameter ;
                :has_value ?value ;
                :has_unit ?unit .
                         
        ?p      :has_parameter ?param ;
                :has_name ?specimen ;
                a :DE.
  
optional{?param :determined_with ?m.		
        ?m      :has_name ?derivation.
    optional{?m :has_input ?e .	
            ?e  :has_name ?input;
            ?e  :has_value ?input_value ;
                :has_unit ?input_unit .}}
}
```