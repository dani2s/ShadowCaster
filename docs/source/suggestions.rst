Suggestions
===========

Installation of OrthoMCL Pipeline
---------------------------------

For the database setup `step <https://github.com/apetkau/orthomcl-pipeline/blob/master/INSTALL.md#step-3-database-setup>`_ , you can previously create the database 
by logging into the MySQL server as root and running:

 ::
			
	mysql> CREATE DATABASE orthomcl;
	mysql> CREATE USER orthomcl@'localhost';
	mysql> set password for orthomcl@localhost = password('yourpassword');
	mysql> GRANT SELECT,INSERT,UPDATE,DELETE,CREATE VIEW,CREATE, INDEX, DROP on orthomcl.* TO orthomcl@localhost;
 
It is important that you create a NEW DATABASE called orthomcl.
