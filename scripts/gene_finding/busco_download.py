import os
import shutil
from os.path import join, exists

from src.logger import print_info
from src.qutils import get_dir_for_download, download_file


def download_db(db_url, clade):
    dirpath = get_dir_for_download('busco', 'busco_db', [clade])
    if not dirpath:
        return None

    db_dirpath = join(dirpath, clade)

    if not exists(db_dirpath):
        downloaded_fpath = join(dirpath, clade + '.tar.gz')
        print_info('  Downloading BUSCO database...')
        download_unpack_compressed_tar(clade + ' database', db_url, downloaded_fpath, db_dirpath)


def download_unpack_compressed_tar(name, download_url, downloaded_fpath, final_dirpath, ext='gz'):
    if download_file(download_url, downloaded_fpath, name, move_file=True):
        unpack_tar(downloaded_fpath, final_dirpath, ext=ext)
        print_info('  Done')
        return True
    return False


def unpack_tar(fpath, dst_dirpath, ext='bz2'):
    import tarfile
    tar = tarfile.open(fpath, "r:" + ext)
    tar.extractall(dst_dirpath)
    tar.close()
    temp_dirpath = join(dst_dirpath, tar.members[0].name)
    from distutils.dir_util import copy_tree
    copy_tree(temp_dirpath, dst_dirpath)
    shutil.rmtree(temp_dirpath)
    os.remove(fpath)
    return True


def main():
    urls = ['https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.2020-03-06.tar.gz',
            'https://busco-data.ezlab.org/v5/data/lineages/fungi_odb10.2021-06-28.tar.gz',
            'https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2020-09-10.tar.gz']
    clades = ['bacteria_odb10', 'fungi_odb10', 'eukaryote_odb10']
    for url, clade in zip(urls, clades):
        download_db(url, clade)
        break


if __name__ == '__main__':
    main()
