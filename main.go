/*
 * Copyright (c) 2019, Aditya Rohan
 * Copyright (c) 2019, Aniket Pandey
 *
 * Submitted to:
 * CS622A: 2019-20 Fall Semester. Course Project
 */

package main

import (
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path"
	"strings"
)

type prefetch struct {
	name          string
	l1d, l2c, llc string
}

// var prefetchers []prefetch

func main() {
	// Initialize info about all prefetchers
	prefetchers := []prefetch{
		{"berti", "berti+", "berti", "next_line_no"},
		{"bingo", "bingo_dpc3", "next_line", "next_line"},
		{"bouquet", "ipcp", "ipcp", "ipcp"},
		{"default", "no", "no", "no"},
		//{"dsp", "dsp", "dsp", "no"},
		{"enhancing", "team_12", "team_12", "team_12"},
		//{"ipcp++", "ipcp++", "cs622_sangam++", "team_12"},
		{"multi-lop", "mlop_dpc3", "next_line", "next_line"},
		{"next-line", "next_line", "next_line", "next_line"},
		{"pangloss", "pangloss", "pangloss", "no"},
		{"sangam", "sangam_dpc3", "sangam_dpc3", "sangam_dpc3"},
		{"sangam++", "cs622_sangam++", "cs622_sangam++", "cs622_sangam++"},
		{"t-skid", "spp", "spp", "no"},
	}

	// Need current directory to build single-core CPUs
	currDir, err := os.Getwd()
	if err != nil {
		log.Fatal(err)
	}
	if len(os.Args) == 1 {
		if err := buildCPUs(currDir, prefetchers); err != nil {
			log.Fatal("Could not build champsim CPUs\n" + err.Error())
		}
	}
	if err := runTraces(currDir, prefetchers); err != nil {
		log.Fatal("Could not run simulation\n" + err.Error())
	}

}

func cleanPrefetchers(prefetchDir string) {
	dir, _ := ioutil.ReadDir(prefetchDir)
	for _, d := range dir {
		os.RemoveAll(path.Join([]string{prefetchDir, d.Name()}...))
	}
}

// Source: https://stackoverflow.com/a/21067803
func copyFileContents(src, dst string) (err error) {
	in, err := os.Open(src)
	if err != nil {
		return
	}
	defer in.Close()
	out, err := os.Create(dst)
	if err != nil {
		return
	}
	defer func() {
		cerr := out.Close()
		if err == nil {
			err = cerr
		}
	}()
	if _, err = io.Copy(out, in); err != nil {
		return
	}
	err = out.Sync()
	return
}

// Build CPUs corresponding to each prefetcher
func buildCPUs(currDir string, prefetchers []prefetch) (err error) {
	prefetchDir := currDir + "/champsim/prefetcher/"
	cleanPrefetchers(prefetchDir)
	for _, prefetcher := range prefetchers {
		// Move prefetching algorithms for L1, L2, LLC to champsim
		fileDir := currDir + "/prefetchers/" + prefetcher.name + "/"
		files, err := ioutil.ReadDir(fileDir)
		if err != nil {
			return err
		}
		for _, file := range files {
			destFile := strings.TrimSuffix(prefetchDir+file.Name(), ".cpp")
			if err := copyFileContents(fileDir+file.Name(),
				destFile); err != nil {
				return err
			}
		}

		coreName := "perceptron-" + prefetcher.l1d + "-" +
			prefetcher.l2c + "-" + prefetcher.llc + "-lru-1core"
		fmt.Println("Building executable: " + coreName)
		cmd := exec.Command("./build_champsim.sh", "perceptron",
			prefetcher.l1d, prefetcher.l2c, prefetcher.llc, "lru", "1")
		cmd.Dir = currDir + "/champsim/"
		if err := cmd.Run(); err != nil {
			return err
		}
		// Remove prefetchers
		for _, file := range files {
			destFile := strings.TrimSuffix(prefetchDir+file.Name(), ".cpp")
			err := os.Remove(destFile)
			if err != nil {
				return err
			}
		}
	}
	return nil
}

// Run Champsim simulator on all traces for all submissions
func runTraces(currDir string, prefetchers []prefetch) (err error) {
	traceDir := currDir + "/champsim/dpc3_traces/"
	files, err := ioutil.ReadDir(traceDir)
	if err != nil {
		log.Fatal(err)
	}
	for _, prefetcher := range prefetchers {
		coreName := "perceptron-" + prefetcher.l1d + "-" +
			prefetcher.l2c + "-" + prefetcher.llc + "-lru-1core"
		for _, traceFile := range files {
			fmt.Println("Running " + coreName + " on " + traceFile.Name())
			cmd := exec.Command("./run_champsim.sh", coreName,
				"1", "10", traceFile.Name())
			cmd.Dir = currDir + "/champsim/"
			if err := cmd.Run(); err != nil {
				return err
			}
		}
	}
	// Move all results to current directory
	resultDir := currDir + "/champsim/results_10M/"
	resultFiles, err := ioutil.ReadDir(resultDir)
	if err != nil {
		return err
	}
	finalDir := currDir + "/results/"
	os.Mkdir(finalDir, os.ModePerm)
	for _, file := range resultFiles {
		if err := copyFileContents(resultDir+file.Name(),
			finalDir+file.Name()); err != nil {
			return err
		}
	}
	fmt.Println("\nAll result files stored in " + finalDir + " directory")
	return nil
}
