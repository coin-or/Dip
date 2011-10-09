#!/usr/local/bin/perl
use strict;
use warnings;

my $MIN_AVAIL = 2;	# how many should i leave for others
my $USER = "magh";
my $EXECUTABLE = "/home/magh/bin/decomp/alps_mad";  #is this visible to all of blaze?
my $SLEEP_RANGE = 60;
my $r = rand($SLEEP_RANGE);
sleep($r);

print("Executing on host:");
system "/bin/hostname";
#system "/usr/local/ilog/ilm/linux/ilmlist";
# use this to see how many i am using
my $my_tokens_in_use = `/usr/local/ilm/linux/ilmlist | grep $USER -c`;
#print("my licences in use = $my_tokens_in_use\n");

my @tokens_available_output = `/usr/local/ilm/linux/ilmlist | grep 'available tokens' `;
chomp $tokens_available_output[2];
print("tokens available = $tokens_available_output[2]\n");
my @tokens_available = ($tokens_available_output[2] =~ m/(\d+)/g);
print("licences available = $tokens_available[0]\n");

if ($tokens_available[0]>$MIN_AVAIL) {
   # This is the successful case
   # start the actual job
   #
   system "$EXECUTABLE @ARGV";
   exit 0;
}
exit 1;

