function first(){
	mhcclassvalue = document.getElementById('mhcclass').value;
	if (mhcclassvalue == 'mhci') {
		document.getElementById('classialleles').style.display = 'block';
		document.getElementById('classiialleles').style.display = 'none';
	}
	else {
		document.getElementById('classialleles').style.display = 'none';
		document.getElementById('classiialleles').style.display = 'block';
	}
}
function second(){
	mhcclassvalue = document.getElementById('mhcclass').value;
	if (mhcclassvalue == 'mhci') {
		document.getElementById('classilength').style.display = 'block';
		document.getElementById('classiilength').style.display = 'none';
	}
	else {
		document.getElementById('classilength').style.display = 'none';
		document.getElementById('classiilength').style.display = 'block';
	}
}
