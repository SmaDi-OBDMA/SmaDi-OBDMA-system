"use strict";

const pages=[
`<div class="w3-padding-16">
Choose a template for the query
</div>

<div class="w3-row-padding" >
<div class="w3-half w3-padding-16">

  <div class="w3-round w3-border w3-hover-shadow" onclick="pageserv.load(1);">
    <header class="w3-container w3-blue">
      <h3>Parameter Query</h3>
    </header>

    <div class="w3-container">
      <p>Querying for smart material parameters</p>
    </div>
  </div>
</div>

<div class="w3-half w3-padding-16">
  
  <div class="w3-round w3-border w3-hover-shadow" onclick="pageserv.load(2);">
    <header class="w3-container w3-blue">
      <h3>Characteristic curve</h3>
    </header>

    <div class="w3-container">
      <p>Querying for characteristic curves of smart materials</p>
    </div>
  </div>

</div>
</div>`,////////////////////////////////////////////////////////////
`<h2> Parameter Query </h2>
  <label for="parameter">Select a parameter from the menu:</label><br>
  <div class="w3-dropdown-click">
    <input type=text class="w3-button w3-blue" id='fclas1' name="isa_result_choice"
        onclick="paramabf.clear();drv.toggleMenu(document.getElementById('dfclas1'));" value="parameter">
    <div class="w3-dropdown-content w3-bar-block w3-card w3-light-grey" id="dfclas1" >
        <input type="text" class="w3-bar-item w3-input" placeholder="searching..." id="sclas1" onkeyup="drv.createMenu('clas1','parameter',this.value);">
        <div id="lisclas1"></div>
    </div>
  </div>
    <br> <br>
    <label for="parameter">Add</label>
    <br>
    <input type="checkbox" id="derivation" name="derivation" onclick="paramabf.tabulate();">
    <label for="derivation"> derivation </label>
    <br>
    <input type="checkbox" id="input" name="input" onclick="paramabf.tabulate();">
    <label for="input"> input </label>
    <br>
    <input type="checkbox" id="description" name="description" onclick="paramabf.tabulate();">
    <label for="description"> description </label>
    <br>
  <br>
  <button class="w3-button w3-round w3-light-blue w3-hover-blue" onclick="paramabf.query();">submit</button> <span id="info"></span>
<br><br>
<div class="" id="filterdiv">
    Filters: <br>
    <ul id="filterlst">
    </ul>
  <button class="w3-button w3-round w3-light-blue w3-hover-blue" onclick="paramabf.remfilter();">remove all filters</button>
</div>
<br>
  <table class="w3-table-all" cellspacing="1" cellpadding="3" bgcolor="#000000" id="partable">
  </table>
  <br>
  <div>SPARQL query:</div>
  <textarea id="spq" name="spq" rows="30" cols="100" disabled="readonly">
  </textarea>
<br>
<br><br>
</div>
<br><br><br><br>
<br><br><br><br>
<br><br><br><br>`,////////////////////////////////////////////////////////////
`<h2> Characteristic Curve Query</h2>
  <label for="param_name">Select a type of characteristic curve from the menu</label><br>
  <div class="w3-dropdown-click">
    <input type=text class="w3-button w3-blue" id='fclas2' name="isa_result_choice"
        onclick="kennabf.visu(0);kennabf.clear();drv.toggleMenu(this.nextSibling.nextSibling);" value="characteristic_curve">
    <div class="w3-dropdown-content w3-bar-block w3-card w3-light-grey" id="dfclas2" >
        <div id="lisclas2"></div>
    </div>
  </div><br><br>
  
  <!--label for="param_name">Choose an element from the menu:</label><br>
  <div class="w3-dropdown-click">
    <input type=text class="w3-button w3-blue" id='fclas3' name="isa_result_choice"
        onclick="drv.toggleMenu(this.nextSibling.nextSibling);" value="elements_of_levels">
    <div class="w3-dropdown-content w3-bar-block w3-card w3-light-grey" id="dfclas3" >
        <div id="lisclas3"></div>
    </div>
  </div><br-->
  <button class="w3-button w3-round w3-light-blue w3-hover-blue" onclick="kennabf.query();kennabf.visu(0);">submit</button> <span id="info"></span>
  
  <div id="visu" class="w3-hide">
  <h3> Visualization </h3>
  <label for="param_name">Specimen:</label>
  <select id="mat" class="w3-button w3-blue" name="mat_choice" onchange="kennabf.clear();kennabf.curvemenu(this.value);">
  </select>

  <br><br>

  <label for="param_name">x-axis:</label> 
  <select id="p1" class="w3-button w3-blue" name="p1_choice" onchange="kennabf.clear();">
      <option value="none">none</option>
  </select> &nbsp

  <label for="param_name">y-axis:</label> 
  <select id="p2" class="w3-button w3-blue" name="p2_choice" onchange="kennabf.clear();">
      <option value="none">none</option>
  </select>
  <br><br>
  <button class="w3-button w3-round w3-light-blue w3-hover-blue" onclick="kennabf.draw();">show</button>
<br>
      <svg id="myPlot" style="width:600px;height:400px"></svg>
      <svg id="legende" style="width:600px;height:200px"></svg>
</div>  
<br>
  <div>SPARQL query:</div>
  <textarea id="spq" name="spq" rows="30" cols="100" disabled="readonly">
  </textarea>
<br>
<br><br><br>`];



const pageserv={
    navmenu:["menu0","menu1","menu2"],
    noinit: function() {},
    initpage: [],
    init() {
        this.initpage=[this.noinit,paramabf.init,kennabf.init];
    },
    load(pagenum) {
        console.log("loading",pagenum);
        for (let i=0;i<this.navmenu.length;i++) {
            let menu=document.getElementById(this.navmenu[i]);
            if (i==pagenum && menu.className.indexOf("w3-blue")<0) {
                menu.className += " w3-blue";
                document.getElementById("content").innerHTML=pages[i];
                this.initpage[i]();
            } else if (i!=pagenum) {
                menu.className = menu.className.replace(" w3-blue", "");
            }
        }
    }
}

const paramabf = {
    filter:{},
    init: function() {
        drv.createMenu("clas1","parameter",0);
    },
    get columns() {
        let cols={
            derivation: document.getElementById("derivation").checked,
            input: document.getElementById("input").checked,
            input_value: document.getElementById("input").checked,
            input_unit: document.getElementById("input").checked,
            description: document.getElementById("description").checked,
        };
        return cols;
    },  
    resp:0,
    addfilter(k,v) {
        this.filter[k]=v;
        this.tabulate();
        let filterlst="";
        for (let key in this.filter) {
            filterlst+="<li>"+key+" = "+this.filter[key]+"</li>\n";
        }
        document.getElementById("filterlst").innerHTML=filterlst;
    },
    remfilter() {
        paramabf.filter={};
        document.getElementById("filterlst").innerHTML="";
        paramabf.tabulate();
    },
    query() {
        let parval=document.getElementById("fclas1").value;
        let qterm=sparql.paramabfq.replace("$isa",parval);
        sparql.ask(qterm,this.tabulate);
    },
    clear() {
        sparql.info("");
        document.getElementById("partable").innerHTML="";
    },
    tabulate(resp) {
        resp=resp?resp:paramabf.resp;
        paramabf.resp=resp;
        let filter= paramabf.filter;
        let cols=paramabf.columns;
        if (!resp || resp.results.bindings.length == 0) {
            sparql.info("No parameter found for this query");
            return;
        };
        let table="  <tr>\n";
        for (let ti of resp.head.vars) {
            if (cols[ti]==undefined || cols[ti]==true) {
                table+="    <th>"+ti+"</th>\n";
            };
        };
        table+="  </tr>\n";
        for (let item of resp.results.bindings) {
            let row="  <tr>\n";
            for (let vi of resp.head.vars) {
                let itemi = item[vi];
                let val=itemi?itemi.value:"";
                if ((filter[vi]==undefined) || filter[vi]==val) {
                    if (cols[vi]==undefined || cols[vi]==true) {
                        row+="    <td onclick='paramabf.addfilter(\""+vi+"\",\""+val+"\");'>"+val+"</td>\n";
                    }
                } else {
                    row=0;
                    break;
                }
            }
            table+=row?row+"  </tr>\n":"";
        };
        document.getElementById("partable").innerHTML=table;
    }
}

const kennabf = {
    init: function() {
        drv.createMenu("clas2","characteristic_curve",0);
    }, 
    query() {
        let kenval=document.getElementById("fclas2").value;
        let qterm=sparql.kennlinieq.replace("$kenn",kenval);
        sparql.ask(qterm,this.visinit);
    },
    alldata:0,
    labels:[0,0],
    legend:{},
    get data() {
        const mat=document.getElementById('mat').value;
        const p1=document.getElementById('p1').value;
        const p2=document.getElementById('p2').value;
        const data={};
        let item0=0;
        for (let item of kennabf.alldata.results.bindings) {
            if (item.specimen.value==mat && item.p1.value==p1 && item.p2.value==p2) {
                let val1=item.w1.value;
                let val2=item.w2.value;
                let vid=item.a.value;
                if(item.p3) {
                    let einheit=(item.e3==undefined || item.e3.value=="dimensionless quantity")?"":     
                            item.e3.value;
                    vid+=item.p3.value+ 
                            "_"+item.w3.value+"_"+einheit;
                }
                if (! data[vid]) {
                    data[vid]=[];
                    if (item.condition) {
                        let einheit=(item.condition_unit==undefined || item.condition_unit.value=="dimensionless quantity")?"":     
                            item.condition_unit.value;
                        kennabf.legend[vid]=""+item.condition.value+ 
                            ": "+item.condition_value.value+" "+einheit;
                    } else if(item.p3) {
                        let einheit=(item.e3==undefined || item.e3.value=="dimensionless quantity")?"":     
                            item.e3.value;
                        kennabf.legend[vid]=""+item.p3.value+ 
                            ": "+item.w3.value+" "+einheit;
                    } else {
                        kennabf.legend[vid]="without boundary condition";
                    }
                }
                if (! isNaN(val1) && ! isNaN(val2)) {
                    data[vid].push([item.w1.value,item.w2.value]);
                    item0=item;
                };
            }
        };
        kennabf.labels[0]=item0["e1"]?item0["e1"].value:"missing unit";
        kennabf.labels[1]=item0["e2"]?item0["e2"].value:"missing unit";
        for (let curve in data) {
            data[curve].sort(function(a,b){return a[0]!=b[0]?a[0]-b[0]:b[1]-a[1]});
        }
        return data;
    },
    proben: {},
    visinit(resp) {
        kennabf.alldata=resp;
        kennabf.proben={};
        if (resp.results.bindings.length==0) {
            sparql.info("No characteristic curves found");
            return;
        };
        let cprobe=resp.results.bindings[0].specimen.value;
        for (let item of resp.results.bindings) {
            if (!kennabf.proben[item.specimen.value]) {
                kennabf.proben[item.specimen.value]={};
            }
            if (!kennabf.proben[item.specimen.value][item.p1.value]) {
                kennabf.proben[item.specimen.value][item.p1.value]={};
            }
            kennabf.proben[item.specimen.value][item.p1.value][item.p2.value]=item.p2.value;
        };
        kennabf.curvemenu(cprobe);
        kennabf.visu(1);
    },  
    visu(show) {
        document.getElementById("visu").className=show?"w3-show":"w3-hide";
    },  
    curvemenu(select_val) {
      let mat_dict=kennabf.proben;
      let mat_select = document.getElementById("mat");
      mat_select.innerHTML="";
      let selected=select_val;
      for(let key in mat_dict) {  
        let option = document.createElement("option");
        option.text = key;
        option.value = key;
        if (selected == key) {
          option.selected = true;
        }
        mat_select.add(option);
      }
      let mat_val = document.getElementById('mat').value;
      kennabf.build_p1_select(mat_val);
      let p1_val = document.getElementById('p1').value;
      kennabf.build_p2_select(p1_val);
      document.getElementById('mat').addEventListener('change', function () {
              select_val = this.value;
              kennabf.build_p1_select(select_val);
              let p1_val = document.getElementById('p1').value
              kennabf.build_p2_select(p1_val);
      });
      document.getElementById('p1').addEventListener('change', function () {
              select_val = this.value;
              kennabf.build_p2_select(select_val);
      });
    },
    clear() {
        sparql.info("");
      d3.selectAll("svg > *").remove();
    },
    draw() {
      d3.selectAll("svg > *").remove();
      const data=kennabf.data;
      let xmin;
      let xmax;
      let ymin;
      let ymax;
      for (let curve in data) {
        let cdata=data[curve];
        if (cdata.length>0) {
            if (xmin==undefined) {
                xmin=cdata[0][0];
                xmax=cdata[0][0];
                ymin=cdata[0][1];
                ymax=cdata[0][1];
            };
            for (let d of cdata) {
                xmin=Math.min(xmin,d[0]);
                xmax=Math.max(xmax,d[0]);
                ymin=Math.min(ymin,d[1]);
                ymax=Math.max(ymax,d[1]);
            };
        };
      };
      let linen=Object.values(data).length;
      const pane=kennabf.getPane(xmin,xmax,ymin,ymax,linen);
      const colors=['Crimson', 'DarkGreen', 'MediumSeaGreen', 'OrangeRed', 'Navy', 'LightSkyBlue', 'RoyalBlue', 'DarkOrange', 'LimeGreen', 'SpringGreen', 'DarkTurquoise', 'LightBlue', 'CadetBlue', 'SeaGree', 'DeepSkyBlue', 'SteelBlue', 'Turquoise', 'Blue', 'LightGreen', 'LightSeaGreen', 'Green', 'DodgerBlue', 'Teal', 'LightSteelBlue', 'CornflowerBlue'];
      let ci=0;
      for (let curve in data) {
        let cdata=data[curve];
        if (cdata.length>0) {
            kennabf.drawOne(pane,cdata,colors[ci%colors.length]);
            kennabf.drawLeg(pane,ci,kennabf.legend[curve],colors[ci%colors.length]);
            ci+=1;
        } else {
            console.log("empty curve",curve);
        };
      };
    },
    getPane(xmin,xmax,ymin,ymax,linen) {
      const label_list=kennabf.labels;
      const ylabel=label_list[1];
      const xlabel=label_list[0];
      const xSize = 600; 
      const ySize = 400;
      const margin = 60;
      const xMax = xSize - margin*2; 
      const yMax = ySize - margin*2; 
      const svg = d3.select("#myPlot")
        .append("svg")
        .attr("width", xSize)
        .attr("height", ySize)
        .append("g")
        .attr("transform","translate(" + margin + "," + margin + ")");
        
      // X Achse
      const x = d3.scaleLinear()
        .domain([xmin, xmax]) 
        .range([0, xMax]);
      svg.append("g")
        .attr("transform", "translate(0," + yMax + ")")
        .call(d3.axisBottom(x));
      svg.append("text")
        .attr("class", "x label")
        .attr("text-anchor", "end")
        .attr("x", xMax)
        .attr("y", yMax)
        .attr("dy", "2.3em")
        .text(xlabel);
        
      // Y Achse
      const y = d3.scaleLinear()
        .domain([ymin, ymax]) 
        .range([ yMax, 0]);

      svg.append("g")
        .call(d3.axisLeft(y))
        
      svg.append("text")
        .attr("class", "y label")
        .attr("text-anchor", "start") //end
        .attr("y", 0)
        .attr("dy", "-0.8em")
        .attr("dx", "-2em")
        .text(ylabel);
    
      console.log("linen",linen);
      const legend = d3.select("#legende")
        .append("svg")
        .attr("width", xSize)
        .attr("height", linen*18)
        .append("g")
        .attr("transform","translate(" + 10 + "," + 10 + ")");
      return {svg:svg,x:x,y:y,legend: legend}
    },
    drawLeg(pane,line,txt,ccolor) {
        const svg=pane.legend;
        svg.append('text')
            .attr("width", 600)
            .attr('x', 0)
            .attr('y', 6+18*line)
            .style("fill", ccolor)
            .text("\u2580 "+txt);
    },
    drawOne(pane,data,ccolor) {
      const svg=pane.svg;
      const x=pane.x;
      const y=pane.y;
      svg.append('g')
        .selectAll("whatever")
        .data(data).enter()
        .append("circle")
        .attr("cx", function (d) { return x(d[0]) } )
        .attr("cy", function (d) { return y(d[1]) } )
        .attr("r", 3)
        .style("fill", ccolor);
      svg.append("path")
        .datum(data)
        .attr("fill", "none")
        .attr("stroke", ccolor)
        .attr("stroke-width", 1.5)
        .attr("d", d3.line()
            .x(function(d) { return x(d[0]) })
            .y(function(d) { return y(d[1]) })
        );
    },
    removeOptions(selectElement) {
         let i, L = selectElement.options.length - 1;
         for(i = L; i >= 0; i--) {
            selectElement.remove(i);
         }
    },
    build_p1_select(select_val) {
        let mat_dict=kennabf.proben;
        let p1_select = document.getElementById("p1");
        let selected=p1_select.value;
        kennabf.removeOptions(p1_select);
        for(let key in mat_dict[select_val]) {
            let option = document.createElement("option");
            option.text = key;
            option.value = key;
            if (selected == key) {
              option.selected = true;
            }
            p1_select.add(option);
        }
      },
    build_p2_select(select_val) {
        let mat_dict=kennabf.proben;
        let mat = document.getElementById('mat').value;
        let p2_select = document.getElementById("p2");
        let selected=p2_select.value;
        kennabf.removeOptions(p2_select);
        for(let key in mat_dict[mat][select_val]) {
            let option = document.createElement("option");
            option.text = mat_dict[mat][select_val][key];
            option.value = mat_dict[mat][select_val][key];
            if (selected == mat_dict[mat][select_val][key]) {
              option.selected = true;
            }
            p2_select.add(option);
        }
      }
}
    
    
pageserv.init();
    
    
    
    
    
    