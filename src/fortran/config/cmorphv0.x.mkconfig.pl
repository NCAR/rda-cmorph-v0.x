#!/usr/bin/perl

# Create configuration files for CMORPH v0.x binary file reads. 
# Command line argument is the date in the format YYYYMMDD

use Time::Local;

my $configdir = "/glade/u/home/tcram/cmorph_v0.x/config";
my ($yyyy, $mm, $dd, $ihr, $epoch);
my $cfile;

$date = shift(@ARGV);
$yyyy = substr($date, 0, 4);
$mm   = substr($date, 4, 2);
$dd   = substr($date, 6, 2);
$cfile = $configdir . "/cmorphv0.x.$date.config";

open(CONFIG,">$cfile");
print CONFIG "LATS -90 90\n";
print CONFIG "LONS 0 360\n";
print CONFIG "HOUR 0 21\n";
print CONFIG "FORM NETCDF\n";
print CONFIG "TIME\n";
for ($ihr=0; $ihr<=21; $ihr+=3) {    
  $epoch = timegm(0,0,$ihr,$dd,$mm-1,$yyyy);
  print CONFIG "$epoch\n";
}
close(CONFIG);