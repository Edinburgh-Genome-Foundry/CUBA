export default {
  list: [
    {
      category: 'Sequence design',
      scenarios: [
        require('./SculptASequence.vue').default,
        require('./DesignOverhangs.vue').default,
        require('./SimulateGGAssemblies').default,
        require('./SimulateMultiMethodAssemblies').default,
        require('./InsertPartsOnBackbones').default,
        require('./SketchConstructs').default,
        require('./ConvertSequenceFiles').default,
        require('./DomesticatePartBatches').default,
        require('./TransferFeatures').default,
        require('./SubstitutePartsOverhangs').default,
        require('./DesignSitesArrays').default
      ]
    },
    {
      category: 'Sequence analysis',
      scenarios: [
        require('./PlotSequenceFeatures.vue').default,
        require('./EvaluateManufacturability.vue').default,
        require('./FindCommonBlocks.vue').default,
        require('./CompareTwoSequences.vue').default,
        require('./RenderSequenticons.vue').default,
        require('./ViewOverhangsCrosstalk.vue').default,
        require('./PredictBadCloneRates.vue').default,
        require('./ExtractFeatures.vue').default
      ]
    },
    {
      category: 'Quality Control',
      scenarios: [
        require('./PredictDigestions').default,
        require('./SelectDigestions').default,
        require('./AnalyzeDigests').default,
        require('./SelectPrimers').default,
        require('./FindSaboteurParts').default,
        require('./DesignPartTestBatches').default
      ]
    },
    {
      category: 'Manufacturing',
      scenarios: [
        require('./CreateAssemblyPicklists.vue').default
        // require('./RearrayPlates.vue').default
      ]
    }
  ]
}
