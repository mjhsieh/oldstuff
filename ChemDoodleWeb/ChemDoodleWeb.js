//
// ChemDoodle Web Components 2.1.0
//
// http://web.chemdoodle.com
//
// Copyright 2009 iChemLabs, LLC.  All rights reserved.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// As a special exception to the GPL, any HTML file which merely makes
// function calls to this code, and for that purpose includes it by
// reference, shall be deemed a separate work for copyright law purposes.
// If you modify this code, you may extend this exception to your version
// of the code, but you are not obligated to do so. If you do not wish to
// do so, delete this exception statement from your version.
//
// As an additional exception to the GPL, you may distribute this
// packed form of the code without the copy of the GPL license normally
// required, provided you include this license notice and a URL through
// which recipients can access the corresponding unpacked source code.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// Please contact iChemLabs <http://www.ichemlabs.com/contact> for
// alternate licensing options.
//

var default_backgroundColor='#FFFFFF';var default_scale=1;var default_rotateAngle=0;var default_bondLength=20;var default_angstromsPerBondLength=1.25;var default_atoms_display=true;var default_atoms_color='#000000';var default_atoms_font_size=12;var default_atoms_font_families=['Helvetica','Arial','Dialog'];var default_atoms_circles=false;var default_atoms_circleDiameter=10;var default_atoms_circleBorderWidth=1;var default_atoms_useJMOLColors=false;var default_bonds_display=true;var default_bonds_color='#000000';var default_bonds_width=1;var default_bonds_saturationWidth=.2;var default_bonds_ends='round';var default_bonds_useJMOLColors=false;var default_bonds_saturationAngle=Math.PI/3;var default_bonds_symmetrical=false;var default_bonds_clearOverlaps=false;var default_bonds_overlapClearWidth=.5;function VisualSpecifications(){this.backgroundColor=default_backgroundColor;this.scale=default_scale;this.rotateAngle=default_rotateAngle;this.bondLength=default_bondLength;this.angstromsPerBondLength=default_angstromsPerBondLength;this.atoms_display=default_atoms_display;this.atoms_color=default_atoms_color;this.atoms_font_size=default_atoms_font_size;this.atoms_font_families=default_atoms_font_families;this.atoms_circles=default_atoms_circles;this.atoms_circleDiameter=default_atoms_circleDiameter;this.atoms_circleBorderWidth=default_atoms_circleBorderWidth;this.atoms_useJMOLColors=default_atoms_useJMOLColors;this.bonds_display=default_bonds_display;this.bonds_color=default_bonds_color;this.bonds_width=default_bonds_width;this.bonds_saturationWidth=default_bonds_saturationWidth;this.bonds_ends=default_bonds_ends;this.bonds_useJMOLColors=default_bonds_useJMOLColors;this.bonds_saturationAngle=default_bonds_saturationAngle;this.bonds_symmetrical=default_bonds_symmetrical;this.bonds_clearOverlaps=default_bonds_clearOverlaps;this.bonds_overlapClearWidth=default_bonds_overlapClearWidth;}
var all=[];var currentPoint=null;var shift=false;var alt=false;function Canvas(){this.lastPoint=null;this.molecule=null;this.emptyMessage=null;this.repaint=function(){var canvas=document.getElementById(this.id);if(canvas.getContext){var ctx=canvas.getContext('2d');ctx.fillStyle=this.specs.backgroundColor;ctx.fillRect(0,0,this.width,this.height);if(this.molecule!=null){ctx.save();ctx.translate(this.width/2,this.height/2);ctx.rotate(this.specs.rotateAngle);ctx.scale(this.specs.scale,this.specs.scale);ctx.translate(-this.width/2,-this.height/2);this.molecule.draw(ctx,this.specs);ctx.restore();}
else
if(this.emptyMessage!=null){ctx.fillStyle='#737683';ctx.textAlign='center';ctx.textBaseline='middle';ctx.font='18px Helvetica';ctx.fillText(this.emptyMessage,this.width/2,this.height/2);}
if(this.drawChildExtras){this.drawChildExtras(ctx);}}}
this.loadMolecule=function(molecule){this.molecule=molecule;this.center();this.molecule.check();this.repaint();}
this.center=function(){var canvas=document.getElementById(this.id);var p=this.molecule.getCenter3D();var center=new Atom('C',this.width/2,this.height/2,0);center.sub3D(p);for(var i=0;i<this.molecule.atoms.length;i++){this.molecule.atoms[i].add3D(center);};var dim=this.molecule.getDimension();this.specs.scale=1;if(dim.x>this.width||dim.y>this.height){this.specs.scale=Math.min(this.width/dim.x,this.height/dim.y)*.9;}}
this.create=function(id,width,height){this.id=id;this.width=width;this.height=height;document.writeln('<canvas class="ChemDoodleWebComponent" id="'+id+'" width="'+width+'" height="'+height+'"></canvas>');this.specs=new VisualSpecifications();}
this.getMolecule=function(){return this.molecule;}
return true;}
function getLocation(element){var objParent=null;var intX=0;var intY=0;do{intX+=element.offsetLeft;intY+=element.offsetTop;objParent=element.offsetParent.tagName;element=element.offsetParent;}
while(objParent!='BODY');return new Point(intX,intY);}
function blockEvent(e){if(e.preventDefault)
e.preventDefault();e.returnValue=false;}
function mouseCoords(ev){if(ev.pageX||ev.pageY){return new Point(ev.pageX,ev.pageY);}
else{return new Point(ev.clientX+document.documentElement.scrollLeft,ev.clientY+document.documentElement.scrollTop);}}
document.onmousedown=function(event){for(var i=0;i<all.length;i++){var canvas=document.getElementById(all[i].id);var p=getLocation(canvas);var insideP=mouseCoords(event);if(insideP.x>p.x&&insideP.y>p.y&&insideP.x<p.x+canvas.width&&insideP.y<p.y+canvas.height){insideP.sub(p);all[i].lastPoint=insideP;if(all[i].mousedown){all[i].mousedown(insideP);}
blockEvent(event);return false;}};}
document.onmouseup=function(event){for(var i=0;i<all.length;i++){all[i].lastPoint=null;if(all[i].mouseup){var canvas=document.getElementById(all[i].id);var p=getLocation(canvas);var insideP=mouseCoords(event);insideP.sub(p);all[i].mouseup(insideP);}};}
document.onmousemove=function(event){currentPoint=mouseCoords(event);for(var i=0;i<all.length;i++){var canvas=document.getElementById(all[i].id);var p=getLocation(canvas);var insideP=mouseCoords(event);if(insideP.x>p.x&&insideP.y>p.y&&insideP.x<p.x+canvas.width&&insideP.y<p.y+canvas.height){insideP.sub(p);if(all[i].lastPoint!=null){if(event.button==2){if(all[i].rightdrag){all[i].rightdrag(insideP);}}
else{if(all[i].drag){all[i].drag(insideP);}}}
else{if(all[i].move){all[i].move(insideP);}}
blockEvent(event);}
else{if(all[i].mouseexit){all[i].mouseexit(insideP);}}};}
document.onkeypress=function(event){for(var i=0;i<all.length;i++){var canvas=document.getElementById(all[i].id);var p=getLocation(canvas);if(currentPoint.x>p.x&&currentPoint.y>p.y&&currentPoint.x<p.x+canvas.width&&currentPoint.y<p.y+canvas.height){if(all[i].keydown&&(event.keyCode==8||event.keyCode==127)){blockEvent(event);return false;}}}}
document.onkeydown=function(event){if(event.keyCode==16){shift=true;}
else
if(event.keyCode==18){alt=true;}
for(var i=0;i<all.length;i++){var canvas=document.getElementById(all[i].id);var p=getLocation(canvas);if(currentPoint.x>p.x&&currentPoint.y>p.y&&currentPoint.x<p.x+canvas.width&&currentPoint.y<p.y+canvas.height){if(all[i].keydown){all[i].keydown(event.keyCode);blockEvent(event);return false;}}}}
document.onkeyup=function(event){if(event.keyCode==16){shift=false;}
else
if(event.keyCode==18){alt=false;}}
function wheel(event){for(var i=0;i<all.length;i++){var canvas=document.getElementById(all[i].id);var p=getLocation(canvas);if(currentPoint.x>p.x&&currentPoint.y>p.y&&currentPoint.x<p.x+canvas.width&&currentPoint.y<p.y+canvas.height){var delta=0;if(!event)
event=window.event;if(event.wheelDelta){delta=event.wheelDelta/120;if(window.opera)
delta=-delta;}
else
if(event.detail){delta=-event.detail/3;}
if(delta){if(all[i].scroll){all[i].scroll(delta);blockEvent(event);break;}}};}}
if(window.addEventListener)
window.addEventListener('DOMMouseScroll',wheel,false);window.onmousewheel=document.onmousewheel=wheel;function DoodleCanvas(id,width,height){all[all.length]=this;this.create(id,width,height);this.specs.atoms_useJMOLColors=true;this.specs.atoms_circleDiameter=7;this.specs.atoms_circleBorderWidth=0;this.isHelp=false;this.helpPos=new Point(this.width-20,20);this.tempAtom=null;this.drawChildExtras=function(ctx){if(this.tempAtom!=null){ctx.strokeStyle='#00FF00';ctx.fillStyle='#00FF00';ctx.lineWidth=1.2;for(var i=0;i<this.molecule.atoms.length;i++){if(this.molecule.atoms[i].isSelected){ctx.beginPath();ctx.moveTo(this.molecule.atoms[i].x,this.molecule.atoms[i].y);ctx.lineTo(this.tempAtom.x,this.tempAtom.y);ctx.stroke();ctx.beginPath();ctx.arc(this.tempAtom.x,this.tempAtom.y,3,0,Math.PI*2,false);ctx.fill();if(this.tempAtom.isOverlap){ctx.strokeStyle='#C10000';ctx.lineWidth=1.2;ctx.beginPath();ctx.arc(this.tempAtom.x,this.tempAtom.y,7,0,Math.PI*2,false);ctx.stroke();}}};}
var radgrad=ctx.createRadialGradient(this.width-20,20,10,this.width-20,20,2);radgrad.addColorStop(0,'#00680F');radgrad.addColorStop(1,'#FFFFFF');ctx.fillStyle=radgrad;ctx.beginPath();ctx.arc(this.helpPos.x,this.helpPos.y,10,0,Math.PI*2,false);ctx.fill();if(this.isHelp){ctx.lineWidth=2;ctx.strokeStyle='black';ctx.stroke();}
ctx.fillStyle=this.isHelp?'red':'black';ctx.textAlign='center';ctx.textBaseline='middle';ctx.font='14px sans-serif';ctx.fillText('?',this.helpPos.x,this.helpPos.y);}
this.drag=function(p){var changed=false;for(var i=0;i<this.molecule.atoms.length;i++){if(this.molecule.atoms[i].isSelected){changed=true;if(p.distance(this.molecule.atoms[i])<7){var x=this.molecule.atoms[i].x;var y=this.molecule.atoms[i].y;var angles=this.molecule.getAngles(this.molecule.atoms[i]);if(angles.length==0){x+=this.specs.bondLength*Math.cos(-Math.PI/6);y+=this.specs.bondLength*Math.sin(-Math.PI/6);}
else
if(angles.length==1){var radian=0;var b=null;for(var j=0;j<this.molecule.bonds.length;j++){if(this.molecule.bonds[j].contains(this.molecule.atoms[i])){b=this.molecule.bonds[j];}};if(b.bondOrder>=3){radian=angles[0]+Math.PI;}
else{var concerned=angles[0]%Math.PI*2;if(isBetween(concerned,0,Math.PI/2)||isBetween(concerned,Math.PI,3*Math.PI/2)){radian=angles[0]+2*Math.PI/3;}
else{radian=angles[0]-2*Math.PI/3;}}
x+=this.specs.bondLength*Math.cos(radian);y-=this.specs.bondLength*Math.sin(radian);}
else{var use=angleBetweenLargest(angles);x+=this.specs.bondLength*Math.cos(use);y-=this.specs.bondLength*Math.sin(use);}
this.tempAtom=new Atom('C',x,y,0);}
else{if(alt&&shift){this.tempAtom=new Atom('C',p.x,p.y,0);}
else{var angle=this.molecule.atoms[i].angle(p);var length=this.molecule.atoms[i].distance(p);if(!shift){length=this.specs.bondLength;}
if(!alt){var increments=Math.floor((angle+Math.PI/12)/(Math.PI/6));angle=increments*Math.PI/6;}
this.tempAtom=new Atom('C',this.molecule.atoms[i].x+length*Math.cos(angle),this.molecule.atoms[i].y-length*Math.sin(angle),0);}}
for(var j=0;j<this.molecule.atoms.length;j++){if(this.molecule.atoms[j].distance(this.tempAtom)<5){this.tempAtom.x=this.molecule.atoms[j].x;this.tempAtom.y=this.molecule.atoms[j].y;this.tempAtom.isOverlap=true;}}
break;}};if(!changed){var dif=new Point(p.x,p.y);dif.sub(this.lastPoint);for(var i=0;i<this.molecule.atoms.length;i++){this.molecule.atoms[i].add(dif);}
for(var i=0;i<this.molecule.rings.length;i++){this.molecule.rings[i].center=this.molecule.rings[i].getCenter();}}
this.lastPoint=p;this.repaint();}
this.mousedown=function(p){if(this.isHelp){window.open('http://web.chemdoodle.com/DoodlerTutorial.html');this.isHelp=false;this.repaint();return;}
for(var i=0;i<this.molecule.atoms.length;i++){if(this.molecule.atoms[i].isHover){this.molecule.atoms[i].isHover=false;this.molecule.atoms[i].isSelected=true;this.drag(p);return;}};for(var i=0;i<this.molecule.bonds.length;i++){if(this.molecule.bonds[i].isHover){this.molecule.bonds[i].isHover=false;this.molecule.bonds[i].bondOrder+=(this.molecule.bonds[i].bondOrder%1)+1;if(this.molecule.bonds[i].bondOrder>3){this.molecule.bonds[i].bondOrder=1;}
this.repaint();return;}};}
this.mouseup=function(p){for(var i=0;i<this.molecule.atoms.length;i++){if(this.tempAtom!=null&&this.molecule.atoms[i].isSelected){if(this.tempAtom.isOverlap){for(var j=0;j<this.molecule.atoms.length;j++){if(this.molecule.atoms[j].distance(this.tempAtom)<5){this.tempAtom=this.molecule.atoms[j];}}}
else{this.molecule.atoms[this.molecule.atoms.length]=this.tempAtom;}
var found=false;for(var j=0;j<this.molecule.bonds.length;j++){if(this.molecule.bonds[j].contains(this.molecule.atoms[i])&&this.molecule.bonds[j].contains(this.tempAtom)){found=true;this.molecule.bonds[j].bondOrder+=(this.molecule.bonds[j].bondOrder%1)+1;if(this.molecule.bonds[j].bondOrder>3){this.molecule.bonds[j].bondOrder=1;}}}
if(!found){this.molecule.bonds[this.molecule.bonds.length]=new Bond(this.molecule.atoms[i],this.tempAtom,1);}
this.molecule.check();}
this.molecule.atoms[i].isSelected=false;};this.tempAtom=null;this.move(p);}
this.move=function(p){var min=Infinity;var hovering=null;for(var i=0;i<this.molecule.atoms.length;i++){this.molecule.atoms[i].isHover=false;var dist=p.distance(this.molecule.atoms[i]);if(dist<this.specs.bondLength&&dist<min){min=dist;hovering=this.molecule.atoms[i];}};for(var i=0;i<this.molecule.bonds.length;i++){this.molecule.bonds[i].isHover=false;var dist=p.distance(this.molecule.bonds[i].getCenter());if(dist<this.specs.bondLength&&dist<min){min=dist;hovering=this.molecule.bonds[i];}};if(hovering!=null){hovering.isHover=true;}
this.isHelp=false;if(p.distance(this.helpPos)<10){this.isHelp=true;}
this.repaint();}
this.keydown=function(key){if(key>=37&&key<=40){var difx=0;var dify=0;if(key==37){difx=-10;}
else
if(key==38){dify=-10;}
else
if(key==39){difx=10;}
else
if(key==40){dify=10;}
for(var i=0;i<this.molecule.atoms.length;i++){this.molecule.atoms[i].x+=difx;this.molecule.atoms[i].y+=dify;}
for(var i=0;i<this.molecule.rings.length;i++){this.molecule.rings[i].center=this.molecule.rings[i].getCenter();}
this.repaint();}
else
if(key==8||key==127){for(var i=0;i<this.molecule.atoms.length;i++){if(this.molecule.atoms[i].isHover){for(var j=0;j<this.molecule.atoms.length;j++){this.molecule.atoms[j].visited=false;}
var connectionsA=[];var connectionsB=[];this.molecule.atoms[i].visited=true;for(var j=0;j<this.molecule.bonds.length;j++){if(this.molecule.bonds[j].contains(this.molecule.atoms[i])){var atoms=[];var bonds=[];var q=new Queue();q.enqueue(this.molecule.bonds[j].getNeighbor(this.molecule.atoms[i]));while(!q.isEmpty()){var a=q.dequeue();if(!a.visited){a.visited=true;atoms[atoms.length]=a;for(var k=0;k<this.molecule.bonds.length;k++){if(this.molecule.bonds[k].contains(a)&&!this.molecule.bonds[k].getNeighbor(a).visited){q.enqueue(this.molecule.bonds[k].getNeighbor(a));bonds[bonds.length]=this.molecule.bonds[k];}}}}
connectionsA[connectionsA.length]=atoms;connectionsB[connectionsB.length]=bonds;}}
var largest=-1;var index=-1;for(var j=0;j<connectionsA.length;j++){if(connectionsA[j].length>largest){index=j;largest=connectionsA[j].length;}}
if(index>-1){this.molecule.atoms=connectionsA[index];this.molecule.bonds=connectionsB[index];this.molecule.check();}
else{var molecule=new Molecule();molecule.atoms[0]=new Atom('C',0,0,0);this.loadMolecule(molecule);}
this.repaint();break;}}}
else{for(var i=0;i<this.molecule.atoms.length;i++){if(this.molecule.atoms[i].isHover){var check=String.fromCharCode(key);var firstMatch=null;var firstAfterMatch=null;var found=false;for(var j=0;j<symbols.length;j++){if(symbols[j]==this.molecule.atoms[i].label){found=true;}
else
if(symbols[j].charAt(0)==check){if(found&&firstAfterMatch==null){firstAfterMatch=symbols[j];}
else
if(firstMatch==null){firstMatch=symbols[j];}}};if(firstAfterMatch!=null){this.molecule.atoms[i].label=firstAfterMatch;}
else
if(firstMatch!=null){this.molecule.atoms[i].label=firstMatch;}
this.molecule.check();this.repaint();break;}}}}
var molecule=new Molecule();molecule.atoms[0]=new Atom('C',0,0,0);this.loadMolecule(molecule);return true;}
DoodleCanvas.prototype=new Canvas();function FileCanvas(id,width,height,action){this.create(id,width,height);form='<br><form name="FileForm" enctype="multipart/form-data" method="POST" action="'+action+'" target="HiddenFileFrame"><input type="file" name="f" /><input type="submit" name="submitbutton" value="Show File" /></form><iframe id="HFF-'+id+'" name="HiddenFileFrame" height="0" width="0" style="display:none;" onLoad="GetMolFromFrame(\'HFF-'+id+'\', '+id+')"></iframe>';document.writeln(form);this.emptyMessage='Click below to load file';this.repaint();return true;}
FileCanvas.prototype=new Canvas();function MolGrabberCanvas(id,width,height,action){this.create(id,width,height);form='<br><form name="MolGrabberForm" method="POST" action="'+action+'" target="HiddenMolGrabberFrame" onSubmit="ValidateMolecule(MolGrabberForm); return false;"><input type="text" name="q" value="" /><input type="submit" name="submitbutton" value="Show Molecule" /></form><iframe id="HMGF-'+id+'" name="HiddenMolGrabberFrame" height="0" width="0" style="display:none;" onLoad="GetMolFromFrame(\'HMGF-'+id+'\', '+id+')"></iframe>';document.writeln(form);this.emptyMessage='Enter search term below';this.repaint();return true;}
MolGrabberCanvas.prototype=new Canvas();function RotatorCanvas(id,width,height,rotate3D){this.create(id,width,height);this.rotate3D=rotate3D;var me=this;var increment=Math.PI/360;this.xIncrement=increment;this.yIncrement=increment;this.zIncrement=increment;this.handle=null;this.timeout=33;this.startRotation=function(){this.stopRotation();this.handle=setInterval(function(){if(me.molecule==null){return;}
if(me.rotate3D){var zrot=$M([[Math.cos(me.zIncrement),-Math.sin(me.zIncrement),0],[Math.sin(me.zIncrement),Math.cos(me.zIncrement),0],[0,0,1]]);var yrot=$M([[Math.cos(-me.yIncrement),0,Math.sin(-me.yIncrement)],[0,1,0],[-Math.sin(-me.yIncrement),0,Math.cos(-me.yIncrement)]]);var xrot=$M([[1,0,0],[0,Math.cos(me.xIncrement),-Math.sin(me.xIncrement)],[0,Math.sin(me.xIncrement),Math.cos(me.xIncrement)]]);for(var i=0;i<me.molecule.atoms.length;i++){var a=me.molecule.atoms[i];var pmvals=$V([a.x-me.width/2,a.y-me.height/2,a.z]);var pm=zrot.x(yrot.x(xrot.x(pmvals)));a.x=pm.e(1)+me.width/2;a.y=pm.e(2)+me.height/2;a.z=pm.e(3);}
for(var i=0;i<me.molecule.rings.length;i++){me.molecule.rings[i].center=me.molecule.rings[i].getCenter();}
if(me.specs.atoms_display&&me.specs.atoms_circles){me.molecule.sortAtomsByZ();}
if(me.specs.bonds_display&&me.specs.bonds_clearOverlaps){me.molecule.sortBondsByZ();}}
else{me.specs.rotateAngle+=me.zIncrement;}
me.repaint();},this.timeout);}
this.stopRotation=function(){if(this.handle!=null){clearInterval(this.handle);this.handle=null;}}
this.isRunning=function(){return this.handle!=null;}
return true;}
RotatorCanvas.prototype=new Canvas();function ViewerCanvas(id,width,height){this.create(id,width,height);return true;}
ViewerCanvas.prototype=new Canvas();var curquat=new Array(4);var lastquat=new Array(4);function TransformCanvas(id,width,height,rotate3D){this.create(id,width,height);all[all.length]=this;this.lastPoint=null;this.rotate3D=rotate3D;this.rightdrag=function(p){var t=new Point(p.x,p.y);t.sub(this.lastPoint);for(var i=0;i<this.molecule.atoms.length;i++){this.molecule.atoms[i].add(t);};this.lastPoint=p;this.repaint();}
this.drag=function(p){if(this.rotate3D==true){trackball(curquat,0,0,0,0);var diameter=Math.max(this.width/4,this.height/4);trackball(lastquat,(this.lastPoint.x-this.width/2)/diameter,-(this.lastPoint.y-this.height/2)/diameter,(p.x-this.width/2)/diameter,-(p.y-this.height/2)/diameter);add_quats(lastquat,curquat,curquat);var mline=new Array(16);build_rotmatrix(mline,curquat);for(var i=0;i<this.molecule.atoms.length;i++){var a=this.molecule.atoms[i];var oldX=-(a.x-this.width/2);var oldY=-(a.y-this.height/2);var oldZ=-a.z;a.x=-(oldX*mline[0]+oldY*mline[1]+oldZ*mline[2])+this.width/2;a.y=-(oldX*mline[4]+oldY*mline[5]+oldZ*mline[6])+this.height/2;a.z=-(oldX*mline[8]+oldY*mline[9]+oldZ*mline[10]);};for(var i=0;i<this.molecule.rings.length;i++){this.molecule.rings[i].center=this.molecule.rings[i].getCenter();}
this.lastPoint=p;if(this.specs.atoms_display&&this.specs.atoms_circles){this.molecule.sortAtomsByZ();}
if(this.specs.bonds_display&&this.specs.bonds_clearOverlaps){this.molecule.sortBondsByZ();}
this.repaint();}
else{var center=new Point(this.width/2,this.height/2);var before=center.angle(this.lastPoint);var after=center.angle(p);this.specs.rotateAngle-=(after-before);this.lastPoint=p;this.repaint();}}
this.scroll=function(delta){this.specs.scale+=delta/100;this.repaint();}
return true;}
TransformCanvas.prototype=new Canvas();var LEEWAY=1.1;function getPointsPerAngstrom(){return default_bondLength/default_angstromsPerBondLength;}
function deduceCovalentBonds(molecule){var pointsPerAngstrom=getPointsPerAngstrom();for(var i=0;i<molecule.atoms.length;i++){for(var j=i+1;j<molecule.atoms.length;j++){var first=molecule.atoms[i];var second=molecule.atoms[j];if(first.distance3D(second)<(covalentRadii[first.label]+covalentRadii[second.label])*pointsPerAngstrom*LEEWAY){molecule.bonds[molecule.bonds.length]=new Bond(first,second,1);}}}}
function removeHydrogens(molecule){var atoms=[];var bonds=[];for(var i=0;i<molecule.bonds.length;i++){if(molecule.bonds[i].a1.label!='H'&&molecule.bonds[i].a2.label!='H'){bonds[bonds.length]=molecule.bonds[i];}}
for(var i=0;i<molecule.atoms.length;i++){if(molecule.atoms[i].label!='H'){atoms[atoms.length]=molecule.atoms[i];}}
molecule.atoms=atoms;molecule.bonds=bonds;}
function Link(data){this.data=data;this.next=null;this.reverse=function(before){if(this.next!=null){this.next.reverse(this);}
this.next=before;}
this.getDataArray=function(array){array[array.length]=data;if(this.next!=null){this.next.getDataArray(array);}}
this.count=function(){if(this.next==null){return 1;}
else{return 1+this.next.count();}}}
function PGraphEdge(i1,i2){if(i1!=null){this.head=new Link(i1);this.head.next=new Link(i2);}
this.getLast=function(){var hold=this.head;while(hold.next!=null){hold=hold.next;}
return hold;}
this.getCopy=function(){var copy=new PGraphEdge();var hold=this.head;var copyHold=new Link(hold.data);copy.head=copyHold;while(hold.next!=null){hold=hold.next;copyHold.next=new Link(hold.data);copyHold=copyHold.next;}
return copy;}
this.merge=function(other){var newPGE=this.getCopy();var same=this.head.data;if(other.head.data!=same&&other.getLast().data!=same){same=this.getLast().data;}
var otherBetweens=other.getCopy();if(newPGE.head.data==same){newPGE.reverse();}
if(other.head.data!=same){otherBetweens.reverse();}
otherBetweens.head=otherBetweens.head.next;newPGE.getLast().next=otherBetweens.head;return newPGE;}
this.connectsTo=function(index){return this.head.data==index||this.getLast().data==index;}
this.isCycle=function(){return this.head.data==this.getLast().data;}
this.size=function(){return this.head.count();}
this.reverse=function(){var last=this.getLast();this.head.reverse(null);this.head=last;}}
function indexOf(array,item){for(var i=0;i<array.length;i++){if(array[i]==item){return i;}};return-1;}
function getRings(molecule){var pGraphEdges=[];var pGraphRings=[];for(var i=0;i<molecule.bonds.length;i++){pGraphEdges[pGraphEdges.length]=new PGraphEdge(indexOf(molecule.atoms,molecule.bonds[i].a1),indexOf(molecule.atoms,molecule.bonds[i].a2));}
while(pGraphEdges.length>0){var counts=new Array(molecule.atoms.length);for(var i=0;i<counts.length;i++){counts[i]=0;};for(var i=0;i<pGraphEdges.length;i++){counts[pGraphEdges[i].head.data]++;counts[pGraphEdges[i].getLast().data]++;};var pick=-1;var min=Infinity;for(var i=0;i<counts.length;i++){if(counts[i]>0&&counts[i]<min){min=counts[i];pick=i;}}
var removing=[];var keep=[];for(var i=0;i<pGraphEdges.length;i++){if(pGraphEdges[i].connectsTo(pick)){removing[removing.length]=pGraphEdges[i];}
else{keep[keep.length]=pGraphEdges[i];}};pGraphEdges=keep;for(var i=0;i<removing.length;i++){for(var j=i+1;j<removing.length;j++){var one=removing[i];var two=removing[j];var newPGE=one.merge(two);var overlap=false;var hold=newPGE.head.next;while(!overlap&&hold!=null){var hold2=hold.next;while(!overlap&&hold2!=null){if(hold2.data==hold.data){overlap=true;}
hold2=hold2.next;}
hold=hold.next;}
if(!overlap){if(newPGE.isCycle()){pGraphRings[pGraphRings.length]=newPGE;}
else{pGraphEdges[pGraphEdges.length]=newPGE;}}}}}
var taken=new Array(pGraphRings.length);var lengths=[];for(var i=0;i<pGraphRings.length;i++){lengths[i]=pGraphRings[i].size();taken[i]=false;};var sssr=molecule.bonds.length-molecule.atoms.length+1;var ringsI=new Array(sssr);for(var i=0;i<sssr;i++){var min=Infinity;var take=-1;for(var j=0;j<lengths.length;j++){if(!taken[j]&&lengths[j]<min){min=lengths[j];take=j;}};ringsI[i]=[];if(take!=-1){pGraphRings[take].head.getDataArray(ringsI[i]);taken[take]=true;}};var rings=new Array(ringsI.length);for(var i=0;i<ringsI.length;i++){var ring=new Ring();for(var j=0;j<ringsI[i].length-1;j++){ring.atoms[j]=molecule.atoms[ringsI[i][j]];}
for(var j=0;j<ring.atoms.length-1;j++){for(var k=0;k<molecule.bonds.length;k++){if(molecule.bonds[k].contains(ring.atoms[j])&&molecule.bonds[k].contains(ring.atoms[j+1])){ring.bonds[ring.bonds.length]=molecule.bonds[k];break;}}}
for(var k=0;k<molecule.bonds.length;k++){if(molecule.bonds[k].contains(ring.atoms[0])&&molecule.bonds[k].contains(ring.atoms[ring.atoms.length-1])){ring.bonds[ring.bonds.length]=molecule.bonds[k];break;}}
rings[i]=ring;}
return rings;}
function copy(molecule){for(var i=0;i<molecule.atoms.length;i++){molecule.atoms[i].metaID=i;}
var newMol=new Molecule();for(var i=0;i<molecule.atoms.length;i++){newMol.atoms[i]=new Atom(molecule.atoms[i].label,molecule.atoms[i].x,molecule.atoms[i].y,molecule.atoms[i].z);}
for(var i=0;i<molecule.bonds.length;i++){newMol.bonds[i]=new Bond(newMol.atoms[molecule.bonds[i].a1.metaID],newMol.atoms[molecule.bonds[i].a2.metaID],molecule.bonds[i].bondOrder);}
return newMol;}
function Point(x,y){this.x=x;this.y=y;this.sub=function(p){this.x-=p.x;this.y-=p.y;}
this.add=function(p){this.x+=p.x;this.y+=p.y;}
this.distance=function(p){return Math.sqrt(Math.pow(p.x-this.x,2)+Math.pow(p.y-this.y,2));}
this.angleForStupidCanvasArcs=function(p){var dx=p.x-this.x;var dy=p.y-this.y;var angle=0;if(dx==0){if(dy==0){angle=0;}else if(dy>0){angle=Math.PI/2;}else{angle=3*Math.PI/2;}}else if(dy==0){if(dx>0){angle=0;}else{angle=Math.PI;}}else{if(dx<0){angle=Math.atan(dy/dx)+Math.PI;}else if(dy<0){angle=Math.atan(dy/dx)+2*Math.PI;}else{angle=Math.atan(dy/dx);}}
while(angle<0){angle+=Math.PI*2;}
angle=angle%(Math.PI*2);return angle;}
this.angle=function(p){var dx=p.x-this.x;var dy=this.y-p.y;var angle=0;if(dx==0){if(dy==0){angle=0;}else if(dy>0){angle=Math.PI/2;}else{angle=3*Math.PI/2;}}else if(dy==0){if(dx>0){angle=0;}else{angle=Math.PI;}}else{if(dx<0){angle=Math.atan(dy/dx)+Math.PI;}else if(dy<0){angle=Math.atan(dy/dx)+2*Math.PI;}else{angle=Math.atan(dy/dx);}}
while(angle<0){angle+=Math.PI*2;}
angle=angle%(Math.PI*2);return angle;}
return true;}
function Ring(){this.atoms=[];this.bonds=[];this.center=null;this.setupBonds=function(){for(var i=0;i<this.bonds.length;i++){this.bonds[i].ring=this;};this.center=this.getCenter();}
this.getCenter=function(){var minX=minY=Infinity;var maxX=maxY=-Infinity;for(var i=0;i<this.atoms.length;i++){minX=Math.min(this.atoms[i].x,minX);minY=Math.min(this.atoms[i].y,minY);maxX=Math.max(this.atoms[i].x,maxX);maxY=Math.max(this.atoms[i].y,maxY);};return new Point((maxX+minX)/2,(maxY+minY)/2);}}
var symbols=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Uub','Uut','Uuq','Uup','Uuh','Uus','Uuo'];var jmolColors=new Array();jmolColors['H']='#ffffff';jmolColors['He']='#d9ffff';jmolColors['Li']='#cc80ff';jmolColors['Be']='#c2ff00';jmolColors['B']='#ffb5b5';jmolColors['C']='#909090';jmolColors['N']='#3050f8';jmolColors['O']='#ff0d0d';jmolColors['F']='#90e050';jmolColors['Ne']='#b3e3f5';jmolColors['Na']='#ab5cf2';jmolColors['Mg']='#8aff00';jmolColors['Al']='#bfa6a6';jmolColors['Si']='#f0c8a0';jmolColors['P']='#ff8000';jmolColors['S']='#ffff30';jmolColors['Cl']='#1ff01f';jmolColors['Ar']='#80d1e3';jmolColors['K']='#8f40d4';jmolColors['Ca']='#3dff00';jmolColors['Sc']='#e6e6e6';jmolColors['Ti']='#bfc2c7';jmolColors['V']='#a6a6ab';jmolColors['Cr']='#8a99c7';jmolColors['Mn']='#9c7ac7';jmolColors['Fe']='#e06633';jmolColors['Co']='#f090a0';jmolColors['Ni']='#50d050';jmolColors['Cu']='#c88033';jmolColors['Zn']='#7d80b0';jmolColors['Ga']='#c28f8f';jmolColors['Ge']='#668f8f';jmolColors['As']='#bd80e3';jmolColors['Se']='#ffa100';jmolColors['Br']='#a62929';jmolColors['Kr']='#5cb8d1';jmolColors['Rb']='#702eb0';jmolColors['Sr']='#ff00';jmolColors['Y']='#94ffff';jmolColors['Zr']='#94e0e0';jmolColors['Nb']='#73c2c9';jmolColors['Mo']='#54b5b5';jmolColors['Tc']='#3b9e9e';jmolColors['Ru']='#248f8f';jmolColors['Rh']='#a7d8c';jmolColors['Pd']='#6985';jmolColors['Ag']='#c0c0c0';jmolColors['Cd']='#ffd98f';jmolColors['In']='#a67573';jmolColors['Sn']='#668080';jmolColors['Sb']='#9e63b5';jmolColors['Te']='#d47a00';jmolColors['I']='#940094';jmolColors['Xe']='#429eb0';jmolColors['Cs']='#57178f';jmolColors['Ba']='#c900';jmolColors['La']='#70d4ff';jmolColors['Ce']='#ffffc7';jmolColors['Pr']='#d9ffc7';jmolColors['Nd']='#c7ffc7';jmolColors['Pm']='#a3ffc7';jmolColors['Sm']='#8fffc7';jmolColors['Eu']='#61ffc7';jmolColors['Gd']='#45ffc7';jmolColors['Tb']='#30ffc7';jmolColors['Dy']='#1fffc7';jmolColors['Ho']='#ff9c';jmolColors['Er']='#e675';jmolColors['Tm']='#d452';jmolColors['Yb']='#bf38';jmolColors['Lu']='#ab24';jmolColors['Hf']='#4dc2ff';jmolColors['Ta']='#4da6ff';jmolColors['W']='#2194d6';jmolColors['Re']='#267dab';jmolColors['Os']='#266696';jmolColors['Ir']='#175487';jmolColors['Pt']='#d0d0e0';jmolColors['Au']='#ffd123';jmolColors['Hg']='#b8b8d0';jmolColors['Tl']='#a6544d';jmolColors['Pb']='#575961';jmolColors['Bi']='#9e4fb5';jmolColors['Po']='#ab5c00';jmolColors['At']='#754f45';jmolColors['Rn']='#428296';jmolColors['Fr']='#420066';jmolColors['Ra']='#7d00';jmolColors['Ac']='#70abfa';jmolColors['Th']='#baff';jmolColors['Pa']='#a1ff';jmolColors['U']='#8fff';jmolColors['Np']='#80ff';jmolColors['Pu']='#6bff';jmolColors['Am']='#545cf2';jmolColors['Cm']='#785ce3';jmolColors['Bk']='#8a4fe3';jmolColors['Cf']='#a136d4';jmolColors['Es']='#b31fd4';jmolColors['Fm']='#b31fba';jmolColors['Md']='#b30da6';jmolColors['No']='#bd0d87';jmolColors['Lr']='#c70066';jmolColors['Rf']='#cc0059';jmolColors['Db']='#d1004f';jmolColors['Sg']='#d90045';jmolColors['Bh']='#e00038';jmolColors['Hs']='#e6002e';jmolColors['Mt']='#eb0026';var covalentRadii=new Array();covalentRadii['H']=0.31;covalentRadii['He']=0.28;covalentRadii['Li']=1.28;covalentRadii['Be']=0.96;covalentRadii['B']=0.84;covalentRadii['C']=0.76;covalentRadii['N']=0.71;covalentRadii['O']=0.66;covalentRadii['F']=0.57;covalentRadii['Ne']=0.58;covalentRadii['Na']=1.66;covalentRadii['Mg']=1.41;covalentRadii['Al']=1.21;covalentRadii['Si']=1.11;covalentRadii['P']=1.07;covalentRadii['S']=1.05;covalentRadii['Cl']=1.02;covalentRadii['Ar']=1.06;covalentRadii['K']=2.03;covalentRadii['Ca']=1.76;covalentRadii['Sc']=1.7;covalentRadii['Ti']=1.6;covalentRadii['V']=1.53;covalentRadii['Cr']=1.39;covalentRadii['Mn']=1.39;covalentRadii['Fe']=1.32;covalentRadii['Co']=1.26;covalentRadii['Ni']=1.24;covalentRadii['Cu']=1.32;covalentRadii['Zn']=1.22;covalentRadii['Ga']=1.22;covalentRadii['Ge']=1.2;covalentRadii['As']=1.19;covalentRadii['Se']=1.2;covalentRadii['Br']=1.2;covalentRadii['Kr']=1.16;covalentRadii['Rb']=2.2;covalentRadii['Sr']=1.95;covalentRadii['Y']=1.9;covalentRadii['Zr']=1.75;covalentRadii['Nb']=1.64;covalentRadii['Mo']=1.54;covalentRadii['Tc']=1.47;covalentRadii['Ru']=1.46;covalentRadii['Rh']=1.42;covalentRadii['Pd']=1.39;covalentRadii['Ag']=1.45;covalentRadii['Cd']=1.44;covalentRadii['In']=1.42;covalentRadii['Sn']=1.39;covalentRadii['Sb']=1.39;covalentRadii['Te']=1.38;covalentRadii['I']=1.39;covalentRadii['Xe']=1.4;covalentRadii['Cs']=2.44;covalentRadii['Ba']=2.15;covalentRadii['La']=2.07;covalentRadii['Ce']=2.04;covalentRadii['Pr']=2.03;covalentRadii['Nd']=2.01;covalentRadii['Pm']=1.99;covalentRadii['Sm']=1.98;covalentRadii['Eu']=1.98;covalentRadii['Gd']=1.96;covalentRadii['Tb']=1.94;covalentRadii['Dy']=1.92;covalentRadii['Ho']=1.92;covalentRadii['Er']=1.89;covalentRadii['Tm']=1.9;covalentRadii['Yb']=1.87;covalentRadii['Lu']=1.87;covalentRadii['Hf']=1.75;covalentRadii['Ta']=1.7;covalentRadii['W']=1.62;covalentRadii['Re']=1.51;covalentRadii['Os']=1.44;covalentRadii['Ir']=1.41;covalentRadii['Pt']=1.36;covalentRadii['Au']=1.36;covalentRadii['Hg']=1.32;covalentRadii['Tl']=1.45;covalentRadii['Pb']=1.46;covalentRadii['Bi']=1.48;covalentRadii['Po']=1.4;covalentRadii['At']=1.5;covalentRadii['Rn']=1.5;covalentRadii['Fr']=2.6;covalentRadii['Ra']=2.21;covalentRadii['Ac']=2.15;covalentRadii['Th']=2.06;covalentRadii['Pa']=2.0;covalentRadii['U']=1.96;covalentRadii['Np']=1.9;covalentRadii['Pu']=1.87;covalentRadii['Am']=1.8;covalentRadii['Cm']=1.69;covalentRadii['Bk']=0;covalentRadii['Cf']=0;covalentRadii['Es']=0;covalentRadii['Fm']=0;covalentRadii['Md']=0;covalentRadii['No']=0;covalentRadii['Lr']=0;covalentRadii['Rf']=0;covalentRadii['Db']=0;covalentRadii['Sg']=0;covalentRadii['Bh']=0;covalentRadii['Hs']=0;covalentRadii['Mt']=0;covalentRadii['Ds']=0;covalentRadii['Rg']=0;covalentRadii['Uub']=0;covalentRadii['Uut']=0;covalentRadii['Uuq']=0;covalentRadii['Uup']=0;covalentRadii['Uuh']=0;covalentRadii['Uus']=0;covalentRadii['Uuo']=0;function Atom(label,x,y,z){this.x=x;this.y=y;this.z=z;this.label=label.replace(/^\s*|\s*$/g,'');this.isLone=false;this.isHover=false;this.isSelected=false;this.isOverlap=false;this.add3D=function(p){this.x+=p.x;this.y+=p.y;this.z+=p.z;}
this.sub3D=function(p){this.x-=p.x;this.y-=p.y;this.z-=p.z;}
this.distance3D=function(p){return Math.sqrt(Math.pow(p.x-this.x,2)+Math.pow(p.y-this.y,2)+Math.pow(p.z-this.z,2));}
this.draw=function(ctx,specs){ctx.fillStyle=specs.atoms_color;if(specs.atoms_useJMOLColors){ctx.fillStyle=jmolColors[this.label];}
if(this.isLone||specs.atoms_circles){ctx.beginPath();ctx.arc(this.x,this.y,specs.atoms_circleDiameter/2,0,Math.PI*2,false);ctx.fill();if(specs.atoms_circleBorderWidth>0){ctx.lineWidth=specs.atoms_circleBorderWidth;ctx.stroke(this.x,this.y,0,Math.PI*2,specs.atoms_circleDiameter/2);}}
else
if(this.isLabelVisible()){var font=specs.atoms_font_size+'px ';for(var i=0;i<specs.atoms_font_families.length;i++){var add=',';if(i==specs.atoms_font_families.length-1){add='';}
font=font+specs.atoms_font_families[i]+add;};ctx.font=font;ctx.textAlign='center';ctx.textBaseline='middle';ctx.fillText(this.label,this.x,this.y);}
if(this.isHover||this.isSelected||this.isOverlap){ctx.strokeStyle=this.isHover?'#885110':'#0060B2';ctx.lineWidth=1.2;ctx.beginPath();ctx.arc(this.x,this.y,7,0,Math.PI*2,false);ctx.stroke();}}
this.isLabelVisible=function(){return this.label!='C';}
return true;}
Atom.prototype=new Point(0,0);function Bond(a1,a2,bondOrder){this.a1=a1;this.a2=a2;this.bondOrder=bondOrder;this.isHover=false;this.ring=null;this.getCenter=function(){return new Point((this.a1.x+this.a2.x)/2,(this.a1.y+this.a2.y)/2);}
this.contains=function(a){return a==this.a1||a==this.a2;}
this.getNeighbor=function(a){return a==this.a1?this.a2:this.a1;}
this.draw=function(ctx,specs){var x1=this.a1.x;var x2=this.a2.x;var y1=this.a1.y;var y2=this.a2.y;var difX=x2-x1;var difY=y2-y1;if(specs.atoms_display&&!specs.atoms_circles&&a1.isLabelVisible()){x1+=difX*.25;y1+=difY*.25;}
if(specs.atoms_display&&!specs.atoms_circles&&a2.isLabelVisible()){x2-=difX*.25;y2-=difY*.25;}
if(specs.bonds_clearOverlaps&&this.a1.distance(this.a2)>specs.bondLength*.8){var xs=x1+difX*.15;var ys=y1+difY*.15;var xf=x2-difX*.15;var yf=y2-difY*.15;ctx.strokeStyle=specs.backgroundColor;ctx.lineWidth=specs.bonds_width+specs.bonds_overlapClearWidth*2;ctx.lineCap='round';ctx.beginPath();ctx.moveTo(xs,ys);ctx.lineTo(xf,yf);ctx.closePath();ctx.stroke();}
ctx.strokeStyle=specs.bonds_color;ctx.lineWidth=specs.bonds_width;ctx.lineCap=specs.bonds_ends;if(specs.bonds_useJMOLColors){var linearGradient=ctx.createLinearGradient(x1,y1,x2,y2);linearGradient.addColorStop(0,jmolColors[this.a1.label]);linearGradient.addColorStop(1,jmolColors[this.a2.label]);ctx.strokeStyle=linearGradient;}
ctx.beginPath();switch(this.bondOrder){case 1:ctx.moveTo(x1,y1);ctx.lineTo(x2,y2);break;case 1.5:ctx.moveTo(x1,y1);ctx.lineTo(x2,y2);break;case 2:if(!specs.bonds_symmetrical&&(this.ring!=null||this.a1.label=='C'&&this.a2.label=='C')){ctx.moveTo(x1,y1);ctx.lineTo(x2,y2);var clip=0;var dist=this.a1.distance(this.a2);var angle=this.a1.angle(this.a2);var perpendicular=angle+Math.PI/2;var useDist=dist*specs.bonds_saturationWidth;var clipAngle=specs.bonds_saturationAngle;if(clipAngle<Math.PI/2){clip=-(useDist/Math.tan(clipAngle));}
if(Math.abs(clip)<dist/2){var xuse1=x1-Math.cos(angle)*clip;var xuse2=x2+Math.cos(angle)*clip;var yuse1=y1+Math.sin(angle)*clip;var yuse2=y2-Math.sin(angle)*clip;var cx1=xuse1-Math.cos(perpendicular)*useDist;var cy1=yuse1+Math.sin(perpendicular)*useDist;var cx2=xuse1+Math.cos(perpendicular)*useDist;var cy2=yuse1-Math.sin(perpendicular)*useDist;var cx3=xuse2-Math.cos(perpendicular)*useDist;var cy3=yuse2+Math.sin(perpendicular)*useDist;var cx4=xuse2+Math.cos(perpendicular)*useDist;var cy4=yuse2-Math.sin(perpendicular)*useDist;var flip=this.ring==null||(this.ring.center.angle(a1)>this.ring.center.angle(a2)&&!(this.ring.center.angle(a1)-this.ring.center.angle(a2)>Math.PI)||(this.ring.center.angle(a1)-this.ring.center.angle(a2)<-Math.PI));if(flip){ctx.moveTo(cx1,cy1);ctx.lineTo(cx3,cy3);}
else{ctx.moveTo(cx2,cy2);ctx.lineTo(cx4,cy4);}}}
else{var useDist=this.a1.distance(this.a2)*specs.bonds_saturationWidth/2;var perpendicular=this.a1.angle(this.a2)+Math.PI/2;var cx1=x1-Math.cos(perpendicular)*useDist;var cy1=y1+Math.sin(perpendicular)*useDist;var cx2=x1+Math.cos(perpendicular)*useDist;var cy2=y1-Math.sin(perpendicular)*useDist;var cx3=x2+Math.cos(perpendicular)*useDist;var cy3=y2-Math.sin(perpendicular)*useDist;var cx4=x2-Math.cos(perpendicular)*useDist;var cy4=y2+Math.sin(perpendicular)*useDist;ctx.moveTo(cx1,cy1);ctx.lineTo(cx4,cy4);ctx.moveTo(cx2,cy2);ctx.lineTo(cx3,cy3);}
break;case 3:var useDist=this.a1.distance(this.a2)*specs.bonds_saturationWidth;var perpendicular=this.a1.angle(this.a2)+Math.PI/2;var cx1=x1-Math.cos(perpendicular)*useDist;var cy1=y1+Math.sin(perpendicular)*useDist;var cx2=x1+Math.cos(perpendicular)*useDist;var cy2=y1-Math.sin(perpendicular)*useDist;var cx3=x2+Math.cos(perpendicular)*useDist;var cy3=y2-Math.sin(perpendicular)*useDist;var cx4=x2-Math.cos(perpendicular)*useDist;var cy4=y2+Math.sin(perpendicular)*useDist;ctx.moveTo(cx1,cy1);ctx.lineTo(cx4,cy4);ctx.moveTo(cx2,cy2);ctx.lineTo(cx3,cy3);ctx.moveTo(x1,y1);ctx.lineTo(x2,y2);break;}
ctx.stroke();if(this.isHover){var angle=this.a1.angleForStupidCanvasArcs(this.a2)+Math.PI/2;ctx.strokeStyle='#885110';ctx.lineWidth=1.2;ctx.beginPath();var angleTo=angle+Math.PI;angleTo=angleTo%(Math.PI*2);ctx.arc(this.a1.x,this.a1.y,6,angle,angleTo,false);ctx.stroke();ctx.beginPath();angle+=Math.PI;angleTo=angle+Math.PI;angleTo=angleTo%(Math.PI*2);ctx.arc(this.a2.x,this.a2.y,7,angle,angleTo,false);ctx.stroke();}}
return true;}
function Molecule(){this.atoms=[];this.bonds=[];this.rings=[];this.findRings=true;this.draw=function(ctx,specs){if(specs.bonds_display==true){for(var i=0;i<this.bonds.length;i++){this.bonds[i].draw(ctx,specs);};}
if(specs.atoms_display==true){for(var i=0;i<this.atoms.length;i++){this.atoms[i].draw(ctx,specs);};}}
this.getCenter3D=function(){if(this.atoms.length==1){return new Atom('C',this.atoms[0].x,this.atoms[0].y,this.atoms[0].z);}
var minX=minY=minZ=Infinity;var maxX=maxY=maxZ=-Infinity;for(var i=0;i<this.atoms.length;i++){minX=Math.min(this.atoms[i].x,minX);minY=Math.min(this.atoms[i].y,minY);minZ=Math.min(this.atoms[i].z,minZ);maxX=Math.max(this.atoms[i].x,maxX);maxY=Math.max(this.atoms[i].y,maxY);maxZ=Math.max(this.atoms[i].z,maxZ);};return new Atom('C',(maxX+minX)/2,(maxY+minY)/2,(maxZ+minZ)/2);}
this.getCenter=function(){if(this.atoms.length==1){return new Point(this.atoms[0].x,this.atoms[0].y);}
var minX=minY=Infinity;var maxX=maxY=-Infinity;for(var i=0;i<this.atoms.length;i++){minX=Math.min(this.atoms[i].x,minX);minY=Math.min(this.atoms[i].y,minY);maxX=Math.max(this.atoms[i].x,maxX);maxY=Math.max(this.atoms[i].y,maxY);};return new Point((maxX+minX)/2,(maxY+minY)/2);}
this.getDimension=function(){if(this.atoms.length==1){return new Point(0,0);}
var minX=minY=Infinity;var maxX=maxY=-Infinity;for(var i=0;i<this.atoms.length;i++){minX=Math.min(this.atoms[i].x,minX);minY=Math.min(this.atoms[i].y,minY);maxX=Math.max(this.atoms[i].x,maxX);maxY=Math.max(this.atoms[i].y,maxY);};return new Point(maxX-minX,maxY-minY);}
this.check=function(){for(var i=0;i<this.atoms.length;i++){this.atoms[i].isLone=false;if(this.atoms[i].label=='C'){var counter=0;for(var j=0;j<this.bonds.length;j++){if(this.bonds[j].a1==this.atoms[i]||this.bonds[j].a2==this.atoms[i]){counter++;}};if(counter==0){this.atoms[i].isLone=true;}}};if(this.findRings){this.rings=getRings(this);for(var i=0;i<this.rings.length;i++){this.rings[i].setupBonds();};}
this.sortAtomsByZ();this.sortBondsByZ();}
this.getAngles=function(a){var angles=[];for(var i=0;i<this.bonds.length;i++){if(this.bonds[i].contains(a)){angles[angles.length]=a.angle(this.bonds[i].getNeighbor(a));}};angles.sort();return angles;}
this.sortAtomsByZ=function(){for(var i=1;i<this.atoms.length;i++){var index=i;while(index>0&&this.atoms[index].z<this.atoms[index-1].z){var hold=this.atoms[index];this.atoms[index]=this.atoms[index-1];this.atoms[index-1]=hold;index--;}}}
this.sortBondsByZ=function(){for(var i=1;i<this.bonds.length;i++){var index=i;while(index>0&&(this.bonds[index].a1.z+this.bonds[index].a2.z)<(this.bonds[index-1].a1.z+this.bonds[index-1].a2.z)){var hold=this.bonds[index];this.bonds[index]=this.bonds[index-1];this.bonds[index-1]=hold;index--;}}}
return true;}
function ValidateMolecule(form){if(form.q.value.match(/^$|^ +$/))
{alert("You must enter a molecule value.");return false;}
form.submit();return true;}
function GetMolFromFrame(frameId,canvas){if(!canvas){return;}
var mol=document.getElementById(frameId).contentDocument.body.innerHTML;if(mol.match('^ChemDoodle Web Components Query Error.')){alert(mol);}
else{var mol=readMOL(mol);removeHydrogens(mol);canvas.loadMolecule(mol);}}
function GetPdbFromFrameDEMO(frameId,canvas){if(!canvas){return;}
var mol=document.getElementById(frameId).contentDocument.body.innerHTML;if(mol.match('^ChemDoodle Web Components Query Error.')){alert(mol);}
else{canvas.loadMolecule(readPDB(mol));}}
function readMOL(content){var molecule=new Molecule();if(content==null||content.length==0){return molecule;}
var currentTagTokens=content.split("\n");var counts=currentTagTokens[3];var numAtoms=parseInt(parseInt(counts.substring(0,3)));var numBonds=parseInt(parseInt(counts.substring(3,6)));for(var i=0;i<numAtoms;i++){var line=currentTagTokens[4+i];molecule.atoms[i]=new Atom(line.substring(31,34),parseFloat(line.substring(0,10))*default_bondLength,-parseFloat(line.substring(10,20))*default_bondLength,parseFloat(line.substring(20,30))*default_bondLength);}
for(var i=0;i<numBonds;i++){var line=currentTagTokens[4+numAtoms+i];var bondOrder=parseInt(line.substring(6,9));if(bondOrder>3){switch(bondOrder){case 4:bondOrder=1.5;break;default:bondOrder=1;break;}}
molecule.bonds[i]=new Bond(molecule.atoms[parseInt(line.substring(0,3))-1],molecule.atoms[parseInt(line.substring(3,6))-1],bondOrder);}
return molecule;}
function writeMOL(molecule){var content='Molecule from ChemDoodle Web Components\n\nhttp://www.ichemlabs.com\n';content=content+fit(molecule.atoms.length.toString(),3)+fit(molecule.bonds.length.toString(),3)+'  0  0  0  0            999 v2000\n';var p=molecule.getCenter();for(var i=0;i<molecule.atoms.length;i++){var a=molecule.atoms[i];content=content+fit(((a.x-p.x)/default_bondLength).toFixed(4),10)+fit((-(a.y-p.y)/default_bondLength).toFixed(4),10)+fit((a.z/default_bondLength).toFixed(4),10)+' '+fit(a.label,3)+' 0  0  0  0  0  0\n';};for(var i=0;i<molecule.bonds.length;i++){var b=molecule.bonds[i];content=content+fit((indexOf(molecule.atoms,b.a1)+1).toString(),3)+fit((indexOf(molecule.atoms,b.a2)+1).toString(),3)+fit(b.bondOrder.toString(),3)+'  0     0  0\n';};content=content+'M  END';return content;}
function fit(data,length){var size=data.length;var padding='';for(var i=0;i<length-size;i++){padding=padding+' ';};return padding+data;}
function readPDB(content){var molecule=new Molecule();if(content==null||content.length==0){return molecule;}
var currentTagTokens=content.split("\n");var pointsPerAngstrom=getPointsPerAngstrom();for(var i=0;i<currentTagTokens.length;i++){var line=currentTagTokens[i];if(line.indexOf("ATOM")==0||line.indexOf("HETATM")==0){molecule.atoms[molecule.atoms.length]=new Atom(line.substring(76,78),parseFloat(line.substring(30,38))*pointsPerAngstrom,-parseFloat(line.substring(46,54))*pointsPerAngstrom,parseFloat(line.substring(38,46)*pointsPerAngstrom));}}
deduceCovalentBonds(molecule);return molecule;}
function angleBetweenLargest(angles){if(angles.length==0){return 0;}
if(angles.length==1){return angles[0]+Math.PI;}
var largest=0;var angle=0;var index=-1;for(var i=0;i<angles.length-1;i++){var dif=angles[i+1]-angles[i];if(dif>largest){largest=dif;angle=(angles[i+1]+angles[i])/2;index=i;}}
var last=angles[0]+Math.PI*2-angles[angles.length-1];if(last>largest){angle=angles[0]-last/2;largest=last;if(angle<0){angle+=Math.PI*2;}
index=angles.length-1;}
return angle;}
function isBetween(x,left,right){return x>=left&&x<=right;}