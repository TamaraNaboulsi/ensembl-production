package VariantsFTP;

################################################################
#
# Connect to different FTP locations that store vcf files;
# based on the division that the organism belongs to
#
################################################################

use Net::FTP;
use Path::Tiny qw(path);
use JSON qw/decode_json/;
use HTTP::Tiny;
use YAML::Tiny;

use Exporter qw(import);
our @EXPORT = qw(connect_to_ftp);

sub connect_to_ftp {
  my ($genome, $release) = @_;
  my $genome_type = $genome->type;
  my $settings_yaml_path = path(__FILE__)->absolute->parent(3)->child('common_files')->child('settings.yml')->stringify;
  my $settings = YAML::Tiny->read($settings_yaml_path)->[0];
  my $ftp_host = $settings->{$genome_type}->{'ftp'}->{'host'};
  my $ftp_directory = get_directory_path($genome, $release, $settings);
  if (!$ftp_directory) {
    print "Something wrong with ftp settings\n";
    return;
  }

  my $ftp_connection = Net::FTP->new($ftp_host, Passive => 1) or die "Can't open $ftp_host\n";
  $ftp_connection->login() or die "Can't log in to the ftp server\n";
  $ftp_connection->cwd($ftp_directory) or return undef;
  $ftp_connection->binary();
  return $ftp_connection;
}

sub get_directory_path {
  my ($genome, $release, $settings) = @_;
  my $ensembl_release = $release->{'ensembl_version'};
  my $ensembl_genomes_release = $release->{'ensembl_genomes_version'};
  my $genome_type = $genome->type;
  my $species = $genome->species();
  my $division = $genome->division;
  my $ftp_pub_dir = $settings->{$genome_type}->{'ftp'}->{'pub_dir'};

  if ($genome_type eq 'vert' || $genome_type eq 'grch37') {
    return "$ftp_pub_dir/release-$ensembl_release/variation/vcf/$species";
  } else {
    # we'll be requesting information from ensembl-genomes ftp
    return "$ftp_pub_dir/release-$ensembl_genomes_release/$division/variation/vcf/$species";
  }
}
