"use strict";

const onto=`[ {
  "@id" : "_:genid1",
  "@type" : [ "http://www.w3.org/2002/07/owl#AllDisjointClasses" ],
  "http://www.w3.org/2002/07/owl#members" : [ {
    "@list" : [ {
      "@id" : "urn:absolute/smadiont#DE"
    }, {
      "@id" : "urn:absolute/smadiont#DE-material"
    }, {
      "@id" : "urn:absolute/smadiont#MSMA"
    }, {
      "@id" : "urn:absolute/smadiont#MSMA-material"
    }, {
      "@id" : "urn:absolute/smadiont#PC"
    }, {
      "@id" : "urn:absolute/smadiont#PC-material"
    }, {
      "@id" : "urn:absolute/smadiont#SMA"
    }, {
      "@id" : "urn:absolute/smadiont#SMA-material"
    } ]
  } ]
}, {
  "@id" : "_:genid10",
  "@type" : [ "http://www.w3.org/2002/07/owl#AllDisjointClasses" ],
  "http://www.w3.org/2002/07/owl#members" : [ {
    "@list" : [ {
      "@id" : "urn:absolute/smadiont#DE"
    }, {
      "@id" : "urn:absolute/smadiont#MSMA"
    }, {
      "@id" : "urn:absolute/smadiont#PC"
    }, {
      "@id" : "urn:absolute/smadiont#SMA"
    } ]
  } ]
}, {
  "@id" : "_:genid15",
  "@type" : [ "http://www.w3.org/2002/07/owl#AllDisjointClasses" ],
  "http://www.w3.org/2002/07/owl#members" : [ {
    "@list" : [ {
      "@id" : "urn:absolute/smadiont#Yeoh_vector"
    }, {
      "@id" : "urn:absolute/smadiont#Yeoh_vector_1"
    }, {
      "@id" : "urn:absolute/smadiont#Yeoh_vector_2"
    }, {
      "@id" : "urn:absolute/smadiont#Yeoh_vector_3"
    }, {
      "@id" : "urn:absolute/smadiont#young_modulus"
    }, {
      "@id" : "urn:absolute/smadiont#young_modulus_neo_Hookean"
    } ]
  } ]
}, {
  "@id" : "_:genid22",
  "@type" : [ "http://www.w3.org/2002/07/owl#AllDisjointClasses" ],
  "http://www.w3.org/2002/07/owl#members" : [ {
    "@list" : [ {
      "@id" : "urn:absolute/smadiont#component_curves"
    }, {
      "@id" : "urn:absolute/smadiont#hysteresis_martensite_content"
    }, {
      "@id" : "urn:absolute/smadiont#magnetic_curve"
    }, {
      "@id" : "urn:absolute/smadiont#material_curves"
    }, {
      "@id" : "urn:absolute/smadiont#stress-strain_curve_elastomer"
    }, {
      "@id" : "urn:absolute/smadiont#stress-strain_curve_martensite"
    } ]
  } ]
}, {
  "@id" : "_:genid29",
  "@type" : [ "http://www.w3.org/2002/07/owl#AllDisjointClasses" ],
  "http://www.w3.org/2002/07/owl#members" : [ {
    "@list" : [ {
      "@id" : "urn:absolute/smadiont#dielectric_strength"
    }, {
      "@id" : "urn:absolute/smadiont#maximum_blocking_stress"
    }, {
      "@id" : "urn:absolute/smadiont#maximum_blocking_stress_hold"
    }, {
      "@id" : "urn:absolute/smadiont#maximum_blocking_stress_load"
    }, {
      "@id" : "urn:absolute/smadiont#relative_permittivity"
    } ]
  } ]
}, {
  "@id" : "urn:absolute/smadiont#",
  "@type" : [ "http://www.w3.org/2002/07/owl#Ontology" ]
}, {
  "@id" : "urn:absolute/smadiont#DE",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#elements_of_materialtypes"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#DE-material",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#DE"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#E-modulus_austenite",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#elastic_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#E-modulus_de-twinned_martensite",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#elastic_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#E-modulus_twinned_martensite",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#elastic_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#MSMA",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#elements_of_materialtypes"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#MSMA-material",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#MSMA"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#MSMA-stick",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#MSMA"
  }, {
    "@id" : "urn:absolute/smadiont#component_level"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#PC",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#elements_of_materialtypes"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#PC-material",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#PC"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#SMA",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#elements_of_materialtypes"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#SMA-material",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#SMA"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#Yeoh_vector",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#elastic_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#Yeoh_vector_1",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#Yeoh_vector"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#Yeoh_vector_2",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#Yeoh_vector"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#Yeoh_vector_3",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#Yeoh_vector"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#austenite_finish_temperature",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#phase_transformation_temperature"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#austenite_start_temperature",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#phase_transformation_temperature"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#characteristic_curve",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ]
}, {
  "@id" : "urn:absolute/smadiont#component_curves",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#characteristic_curve"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#component_level",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#elements_of_levels"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#continuum_level",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#elements_of_levels"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#corresponds_with",
  "@type" : [ "http://www.w3.org/2002/07/owl#ObjectProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topObjectProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#curve_magnetic_stress",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_curves"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#curve_param1",
  "@type" : [ "http://www.w3.org/2002/07/owl#ObjectProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topObjectProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#curve_param2",
  "@type" : [ "http://www.w3.org/2002/07/owl#ObjectProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topObjectProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#determined_with",
  "@type" : [ "http://www.w3.org/2002/07/owl#ObjectProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topObjectProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#dielectric_strength",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#elastic_parameters",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#electric_field",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#state_dependent_quantity"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#electric_permittivity_parameters",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#electrostatic_pressure_model",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#model"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#elements_of_levels",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ]
}, {
  "@id" : "urn:absolute/smadiont#elements_of_materialtypes",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ]
}, {
  "@id" : "urn:absolute/smadiont#flux_density",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#state_dependent_quantity"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#hardening_parameter",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_boundary_condition",
  "@type" : [ "http://www.w3.org/2002/07/owl#ObjectProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topObjectProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_column",
  "@type" : [ "http://www.w3.org/2002/07/owl#DatatypeProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topDataProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_curve",
  "@type" : [ "http://www.w3.org/2002/07/owl#ObjectProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topObjectProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_description",
  "@type" : [ "http://www.w3.org/2002/07/owl#DatatypeProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topDataProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_equation",
  "@type" : [ "http://www.w3.org/2002/07/owl#DatatypeProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topDataProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_input",
  "@type" : [ "http://www.w3.org/2002/07/owl#ObjectProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topObjectProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_input_characteristic_curve",
  "@type" : [ "http://www.w3.org/2002/07/owl#ObjectProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "urn:absolute/smadiont#has_input"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_input_parameter",
  "@type" : [ "http://www.w3.org/2002/07/owl#ObjectProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "urn:absolute/smadiont#has_input"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_name",
  "@type" : [ "http://www.w3.org/2002/07/owl#DatatypeProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topDataProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_output",
  "@type" : [ "http://www.w3.org/2002/07/owl#ObjectProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topObjectProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_parameter",
  "@type" : [ "http://www.w3.org/2002/07/owl#ObjectProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topObjectProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_row",
  "@type" : [ "http://www.w3.org/2002/07/owl#DatatypeProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topDataProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_symbol",
  "@type" : [ "http://www.w3.org/2002/07/owl#DatatypeProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topDataProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_tuple",
  "@type" : [ "http://www.w3.org/2002/07/owl#ObjectProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topObjectProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_unit",
  "@type" : [ "http://www.w3.org/2002/07/owl#DatatypeProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topDataProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#has_value",
  "@type" : [ "http://www.w3.org/2002/07/owl#DatatypeProperty" ],
  "http://www.w3.org/2000/01/rdf-schema#subPropertyOf" : [ {
    "@id" : "http://www.w3.org/2002/07/owl#topDataProperty"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#hysteresis_martensite_content",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#component_curves"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#identification_parameters",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#initial_strain",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#identification_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#initial_stress",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#identification_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#inverse_permittivity_at_constant_mechanic_strain",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#electric_permittivity_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#inverse_permittivity_at_constant_mechanic_stress",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#electric_permittivity_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#lattice_constant",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#lattice_constant_a",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#lattice_constant"
  } ],
  "http://www.w3.org/2002/07/owl#disjointWith" : [ {
    "@id" : "urn:absolute/smadiont#lattice_constant_c"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#lattice_constant_c",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#lattice_constant"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#mag_energy_density",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#state_dependent_quantity"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#magnetic_curve",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_curves"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#magnetic_curve_easy_axis",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#magnetic_curve"
  } ],
  "http://www.w3.org/2002/07/owl#disjointWith" : [ {
    "@id" : "urn:absolute/smadiont#magnetic_curve_hard_axis"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#magnetic_curve_hard_axis",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#magnetic_curve"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#magnetic_field_strength",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#state_dependent_quantity"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#magnetic_polarization",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#state_dependent_quantity"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#magnetic_stress",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#stress_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#magnetization",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#state_dependent_quantity"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#martensite_content",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#state_dependent_quantity"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#martensite_finish_temperature",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#phase_transformation_temperature"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#martensite_start_temperature",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#phase_transformation_temperature"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#material",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#micro_level"
  } ],
  "http://www.w3.org/2002/07/owl#disjointWith" : [ {
    "@id" : "urn:absolute/smadiont#sample"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#material_curves",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#characteristic_curve"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#material_parameter",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#max_magnetic_stress",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#maximum_blocking_stress",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#maximum_blocking_stress_hold",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#maximum_blocking_stress"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#maximum_blocking_stress_load",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#maximum_blocking_stress"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#maximum_strain_mag",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#strain_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#measurement",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ]
}, {
  "@id" : "urn:absolute/smadiont#mechanic_stress",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#stress_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#mechanical_compliance",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#mechanical_compliance_at_constant_electric_field",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#mechanical_compliance"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#mechanical_compliance_at_constant_flux_density",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#mechanical_compliance"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#mechanical_stiffness",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#mechanical_stiffness_at_constant_electric_field",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#mechanical_stiffness"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#mechanical_stiffness_at_constant_flux_density",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#mechanical_stiffness"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#micro_level",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#elements_of_levels"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#model",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ]
}, {
  "@id" : "urn:absolute/smadiont#parameter",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#physical_quantity"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#permittivity_at_constant_mechanic_strain",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#electric_permittivity_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#permittivity_at_constant_mechanic_stress",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#electric_permittivity_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#phase_transformation_temperature",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#physical_constant",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#physical_quantity",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ]
}, {
  "@id" : "urn:absolute/smadiont#piezoelectric_charge_constant",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#piezoelectric_coupling_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#piezoelectric_coupling",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#piezoelectric_coupling_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#piezoelectric_coupling_parameters",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#piezoelectric_strain_constant",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#piezoelectric_coupling_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#piezoelectric_stress_constant",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#piezoelectric_coupling_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#plateau_finish_strain",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#strain_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#plateau_start_strain",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#strain_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#relative_permeability",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#state_dependent_quantity"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#relative_permittivity",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#sample",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#micro_level"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#state_dependent_quantity",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#physical_quantity"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#strain",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#state_dependent_quantity"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#strain_parameter",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#stress-strain_curve_elastomer",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_curves"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#stress-strain_curve_martensite",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_curves"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#stress_parameter",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#state_dependent_quantity"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#temperature",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#state_dependent_quantity"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#tensile_test",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#measurement"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#twinning_stress",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#material_parameter"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#vacuum_permeability",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#physical_constant"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#vacuum_permittivity",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#physical_constant"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#wire_actuator",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#SMA"
  }, {
    "@id" : "urn:absolute/smadiont#component_level"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#yeoh_minR",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#model"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#young_modulus",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#elastic_parameters"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#young_modulus_Neo_Hookean_initial_gradient",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#model"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#young_modulus_Neo_Hookean_minR",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#model"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#young_modulus_initial_gradient",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#model"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#young_modulus_minR",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#model"
  } ]
}, {
  "@id" : "urn:absolute/smadiont#young_modulus_neo_Hookean",
  "@type" : [ "http://www.w3.org/2002/07/owl#Class" ],
  "http://www.w3.org/2000/01/rdf-schema#subClassOf" : [ {
    "@id" : "urn:absolute/smadiont#elastic_parameters"
  } ]
} ]`;