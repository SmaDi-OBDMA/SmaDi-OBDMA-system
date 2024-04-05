"use strict";



// model
const drm = {
    clont:0,
    clidx:0,
    get onto() {
        if (! this.clon) this.initclon();
        return clon;
    },
    initclon() {
        let ont0=onto.replace(/\@/mg,'').replace(/["]([a-zA-Z0-9-]|:|\.|\/)*#/mig, '"');
        ont0=JSON.parse(ont0);
        this.clont={'clas':{},'prop':{}};
        this.clidx={'clas':[],'prop':[]};
        let subnm={'clas':'subClassOf','prop':'subPropertyOf'};
        for (let obj of ont0) {
            let ctyp=obj.type.indexOf("Class")>=0?'clas':
                (obj.type.indexOf("ObjectProperty")>=0?'prop':0);
            if (ctyp) {
                let opart=this.clont[ctyp];
                if (! opart[obj.id]) {
                        opart[obj.id]=[];
                };
                if (obj[subnm[ctyp]]) {
                    for (let subs of obj[subnm[ctyp]]) {
                        if (! opart[subs.id]) {
                            opart[subs.id]=[];
                        }
                        opart[subs.id].push(obj.id);
                    }
                };
            };
        };
        console.log("ontoclon initialised",
            Object.keys(this.clont['clas']).length, "classes,",
            Object.keys(this.clont['prop']).length, "properties");
        for (let ctyp in this.clont) {
            let opart=this.clont[ctyp];
            for (let obj in opart) {
                opart[obj].sort(); 
            };
            this.clidx[ctyp]=Object.keys(this.clont[ctyp]).sort();
        };
    },
    subGen(typ,objid,filter) {
        return this.clont[typ][objid];
    }
}

//view
const drv = {
    caret:' <i class="fa fa-caret-down"></i>',
    submenus(outer,item,filter,subs,d=0,accu='') {
        let cbname='"drv.selectName(\''+outer+'\',\''+item+'\');"';
        let cbcaret='"drv.toggleMenu(this.parentNode.nextSibling);"';
        let cbnbutton='<span onclick='+cbname+'>';
        let cbcbutton='<span onclick='+cbcaret+'>';
        if (subs[item].length>0 && d<5) {
            accu+='<div class="w3-dropdown-click">';
            accu+='<div class="w3-bar w3-button">';
            accu+=cbnbutton+item+'</span>';//'\');">'+
            accu+=cbcbutton+this.caret+'</span>';
            accu+='</div>';
            accu+='<div class="w3-dropdown-content w3-bar-block w3-card w3-light-grey w3-margin-left">';
            let foundsome=false;
            for (let si of subs[item]) {
                let retval=this.submenus(outer,si,filter,subs,d+1,accu);
                if (retval) {
                    accu=retval;
                    foundsome=true;
                };
            };
            if (foundsome) {
                return accu+'</div></div>';
            } else {
                return 0;
            };
        } else {
            if (!filter || item.indexOf(filter) >= 0) {
                return accu+'<span class="w3-bar-item w3-button"' + 
                    'onclick='+cbname+'>'+item+'</span>';
            } else {
                return 0;
            };
        };
    },
    selectName(outer,val) {
        document.getElementById('f'+outer).value=val;
        let mainm=document.getElementById('f'+outer);
        drv.toggleMenu(mainm.nextSibling.nextSibling);
    },
    remSelMenu(outer) { 
        let sbox=document.getElementById("df"+outer);
        sbox.innerHTML="";
        sbox.style.border="0px";
    },
    createMenu(outer,str,filter=0) {
        let cbfname='drv.toggleMenu(this);';
        let sbox=document.getElementById("df"+outer);
        let menucnt='';
        let ctyp=outer.indexOf('clas')>=0?'clas':'prop';
        str=str?str:'';
        let itemlist=drm.subGen(ctyp,str,filter);
        let addit=0;
        let subs=outer.indexOf('clas')>=0?drm.clont['clas']:drm.clont['prop'];
        for (let item of itemlist) {
            let subits=subs[item];
            addit=this.submenus(outer,item,filter,subs);
            menucnt+=addit?addit:"";
        };
        document.getElementById("lis"+outer).innerHTML=menucnt;
    },
    toggleMenu(menu) {
      if (menu.className.indexOf("w3-show") == -1) {
        menu.className += " w3-show";
      } else {
        menu.className = menu.className.replace(" w3-show", "");
      }
    }
}

const sparql={
    url:"http://localhost:8080/sparql",
    prefix:`PREFIX : <urn:absolute/smadiont#> 
        PREFIX owl: <http://www.w3.org/2002/07/owl#>`,
    info: function(txt) {
        let info=document.getElementById("info");
        info.innerHTML=txt;
    },
    ask(query,cbf) {
        sparql.info("Query submitted, please wait.");
        let spq=document.getElementById("spq");
        spq.value=query;
        const xhttp = new XMLHttpRequest();
        xhttp.onreadystatechange = function() {
            if (this.readyState == 4) {
                if (this.status == 200) {
                    sparql.info("");
                    let resp=this.responseText;
                    //console.log("resp",resp);
                    resp=resp.replaceAll("\n\"","\"");
                    try {
                      cbf(JSON.parse(resp));
                    } catch (e) {
                      if (e instanceof SyntaxError){ 
                        console.log(e.message);
                        console.log(e.name); // "SyntaxError"
                        console.log(e.stack); // Stack of the error
                        sparql.info("SPARQL server delivered invalid response.");
                        }
                    }
                } else {
                    console.log("Network returned",this.readyState,this.status,this.statusText);
                    sparql.info("SPARQL server is unreachable.");
                }
            }
        };
        xhttp.open("POST", this.url);
        xhttp.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
        xhttp.send("query=" + encodeURIComponent(this.prefix+query) + "&format=json");
    },
    matq:`
        Select ?Material
        Where {?m a :Material ;           
             :hat_Name ?Material .}`,
    paramq:`
        Select ?parameter
        Where {?m a :parameter ;
                :has_name ?parameter .}`,
    elemq:`
        Select ?specimen
        Where {?m a :elements_of_levels ;
                :has_name ?specimen .}`,
    kennq:`
        Select ?curve
        Where {?m a :characteristic_curve ;
                :has_name ?curve .}`,
    paramabfq:`
    SELECT ?material_class ?specimen  ?parameter  ?value ?unit ?derivation 
    ?input ?input_value ?input_unit ?metadata
WHERE { ?param  a  :$isa ; 	
            :has_name ?parameter ;
            :has_value ?value ;
            :has_unit ?unit .
                         
     ?p     :has_parameter ?param ;
           :has_name ?specimen . 
    values (?type ?material_class) {
            (:MSMA "MSMA")
            (:DE "DE")
            (:SMA "SMA")
            (:PC "PC")}
            ?p a ?type.  
     optional{?param :determined_with ?m.		
             ?m     :has_name ?derivation.
         optional{?m  :has_input ?e .	
                  ?e   :has_name ?input.
         optional{?e   :has_value ?input_value ;
                       :has_unit ?input_unit .}}
         optional{?m :has_description ?metadata .}}
}`,
    kennlinieq:`
SELECT ?a ?p1 ?w1 ?e1 ?p2  ?w2 ?e2 ?p3 ?w3 ?e3 ?specimen ?condition ?condition_unit ?condition_value
WHERE { ?a      a :$kenn.	 
  		?a      :has_tuple ?x .		
		?x      :curve_param1 ?param1.	
		?x      :curve_param2 ?param2.
		?param1 :has_name ?p1;
		        :has_value ?w1;
		        :has_unit ?e1 .
		?param2 :has_name ?p2;
		        :has_value ?w2 ;
		        :has_unit ?e2 .
		
  		?mat    :has_parameter ?param1 .
		?mat    :has_name ?specimen .
		optional{ ?a :has_boundary_condition ?cond .
                    ?cond :has_name ?condition;
                        :has_unit ?condition_unit;
                        :has_value ?condition_value.}
}`
}

drm.initclon();


