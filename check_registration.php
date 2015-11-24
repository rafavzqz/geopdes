

<!DOCTYPE HTML>
<!--
	Eleganta by TEMPLATED
    templated.co @templatedco
    Released for free under the Creative Commons Attribution 3.0 license (templated.co/license)
-->
<html>
	<head>
		<title>Lions-Magenes Days</title>
		<meta http-equiv="content-type" content="text/html; charset=utf-8" />
		<meta name="description" content="" />
		<meta name="keywords" content="" />
		<meta name="viewport" content="width=1040" />
		<link href="http://fonts.googleapis.com/css?family=Open+Sans:400,300,600,700,800" rel="stylesheet" type="text/css" />
		<!--[if lte IE 8]><script src="js/html5shiv.js"></script><![endif]-->
		<script src="http://ajax.googleapis.com/ajax/libs/jquery/1.11.0/jquery.min.js"></script>
		<script src="js/skel.min.js"></script>
		<script src="js/skel-panels.min.js"></script>
		<script src="js/init.js"></script>
		<noscript>
			<link rel="stylesheet" href="css/skel-noscript.css" />
			<link rel="stylesheet" href="css/style.css" />
			<link rel="stylesheet" href="css/style-desktop.css" />
		</noscript>
		<!--[if lte IE 8]><link rel="stylesheet" href="css/ie/v8.css" /><![endif]-->
		<!--[if lte IE 9]><link rel="stylesheet" href="css/ie/v9.css" /><![endif]-->
	</head>
	<body class="homepage">

		<div id="wrapper">

			<!-- Header -->
			<div id="header">
				<div class="container">
					
					<!-- Logo -->
					<div id="logo">
						<h1><a href="#">Lions-Magenes Days</a></h1>
						<span>Pavia, April 13-14, 2015</span>
					</div>
					
				</div>
			</div>
			<!-- Header -->
			
			<div class="container">
				<!-- Nav -->
				<?php
				
				$name=$_POST['name'];
				$position=$_POST['position'];
				$institution=$_POST['institution'];
				$email=$_POST['email'];
				$hotel=$_POST['hotel'];
				
				include('menu.php');
				
				?>
				<!-- /Nav -->
			</div>
			
			<div class="container">
				<div id="main">
					<!-- Banner -->
					<div id="banner">
						<a href="#" class="image full"><img src="images/pavia033.jpg" alt=""></a>
					</div>
					<!-- /Banner -->
					<div class="row">
					<div class="3u">
					<img src="images/Lions.jpg" width="200" height="320" alt="" border="0">
					</div>	
					<div class="6u" align='justify'>
						<section>
							
							
							
					Registration is free but mandatory. 
					
					<hr>
					<h3>Register online</h3>
					
					<br>
					
					<br>
					
					

<table>

<tr>
       <td>Name & Surname:</td>
       <td> <?php echo $name; ?></td>
</tr>
<tr>
       <td>Position</td>
       <td> <?php echo $position; ?></td>
</tr>
<tr>
       <td>Institution</td>
       <td> <?php echo $institution; ?> </td>
</tr>
<tr>
       <td>E-mail address</td>
       <td><?php echo $email; ?></td>
</tr>
<tr>
    </tr>
</table>
<p> 
 <?php if(isset($_POST['hotel'])) 
  	   		{
			switch($hotel){
			case "1":
        echo "You wish to reserve a room at Collegio Borromeo";
        break;
    case "2":
        echo "You wish to reserve a room at Hotel Excelsior";
        break;  
    default:
        echo "You do not wish to reserve a room";
			}}
   
				 else
  { $hotel=0;
   echo 'You do not wish to receive information on accomodation';}
?>

<p> Are these information correct?</p>

<table>

<tr>
       <td><form action="register.php" method="get">
<input type="hidden" name="hotel" value="<?php echo $hotel; ?>">
<input type="hidden" name="name" value="<?php echo $name; ?>">
<input type="hidden" name="position" value="<?php echo $position; ?>">
<input type="hidden" name="institution" value="<?php echo $institution; ?>">
<input type="hidden" name="email" value="<?php echo $email; ?>">
<input type="submit" value="Yes, register">&nbsp;  &nbsp; 
</form> </td>
       <td><form action="registration.php" method="get">
<input type="hidden" name="hotel" value="<?php echo $hotel; ?>">
<input type="hidden" name="name" value="<?php echo $name; ?>">
<input type="hidden" name="position" value="<?php echo $position; ?>">
<input type="hidden" name="institution" value="<?php echo $institution; ?>">
<input type="hidden" name="email" value="<?php echo $email; ?>">
<input type="submit" value="No, back to registration form">
</form></td>
</tr>
</table>



					
				
						
							</section>
				</div>
					
					
					
					<div class="3u">
					<img align='right' src="images/Magenes.jpg" width="200" height="320" alt="" border="0">
					</div>
			
						
					</div> 
					
				</div>
			</div>
	<?php
	
	 include('footer.php'); 
	
	?>
		</div>

		<!-- Copyright -->
		<div id="copyright">
			<div class="container">
				Design: <a href="http://templated.co">TEMPLATED</a> Images: <a href="http://unsplash.com">Unsplash</a> (<a href="http://unsplash.com/cc0">CC0</a>)
			</div>
		</div>

	</body>
</html>





?>