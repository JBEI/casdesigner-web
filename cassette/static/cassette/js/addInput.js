var counter = 1;
var limit = 10;
function addInput(divName){
     if (counter == limit)  {
          alert("You have reached the limit of adding " + counter + " genes");
     }
     else {
          var newdiv = document.createElement('div');
          newdiv.innerHTML = "Gene " + (counter + 1) + "<input type='text' id='geneList1[]' name='geneList1[]'> <br>" +  
          					 "Sequence" + (counter + 1) + "<input type='text' id='seqList1[]' name='seqList1[]''> <br>";
          document.getElementById(divName).appendChild(newdiv);
          counter++;
     }
}